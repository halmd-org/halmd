/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/range/iterator_range.hpp>
#include <exception>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/utility/luabind.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace boost::fusion;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu
{

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
void neighbour<dimension, float_type>::options(po::options_description& desc)
{
    desc.add_options()
        ("cell-occupancy", po::value<float>()->default_value(0.4),
         "desired average cell occupancy")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    using namespace luabind;
    register_any_converter<float>();
}


/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void neighbour<dimension, float_type>::depends()
{
    modules::depends<_Self, particle_type>::required();
    modules::depends<_Self, box_type>::required();
    modules::depends<_Self, force_type>::required();
}

template <int dimension, typename float_type>
neighbour<dimension, float_type>::neighbour(modules::factory& factory, po::variables_map const& vm)
  : _Base(factory, vm)
  // dependency injection
  , particle(modules::fetch<particle_type>(factory, vm))
  , force(modules::fetch<force_type>(factory, vm))
  , box(modules::fetch<box_type>(factory, vm))
  // select thread-dependent reduction kernel
  , dim_reduce(64, (64 << DEVICE_SCALE))
  , displacement_impl(get_displacement_impl(dim_reduce.threads_per_block()))
  // allocate parameters
  , r_skin_(vm["skin"].as<float>())
  , rr_skin_half_(pow(r_skin_ / 2, 2))
  , rr_cut_skin_(particle->ntype, particle->ntype)
  , g_rr_cut_skin_(rr_cut_skin_.data().size())
  , nu_cell_(vm["cell-occupancy"].as<float>())
  , g_r0_(particle->nbox)
  , g_rr_(dim_reduce.blocks_per_grid())
  , h_rr_(g_rr_.size())
  , sort_(particle->nbox, particle->dim.threads_per_block())
{
    /*@{ FIXME remove pre-Lua hack */
    shared_ptr<profiler_type> profiler(modules::fetch<profiler_type>(factory, vm));
    register_runtimes(*profiler);
    /*@}*/

    matrix_type r_cut = force->cutoff();
    typename matrix_type::value_type r_cut_max = 0;
    for (size_t i = 0; i < particle->ntype; ++i) {
        for (size_t j = i; j < particle->ntype; ++j) {
            rr_cut_skin_(i, j) = std::pow(r_cut(i, j) + r_skin_, 2);
            r_cut_max = max(r_cut(i, j), r_cut_max);
        }
    }
    // find an optimal(?) cell size
    // ideally, we would like to have warp_size placeholders per cell
    //
    // definitions:
    // a) n_cell = L_box / cell_length   (for each dimension)
    // b) #cells = prod(n_cell)
    // c) #particles = #cells × cell_size × ν_eff
    //
    // constraints:
    // 1) cell_length > r_c + r_skin  (potential cutoff + neighbour list skin)
    // 2) cell_size is a (small) multiple of warp_size
    // 3) ν_eff ≤ ν  (actual vs. desired occupancy)
    //
    // upper bound on ncell_ provided by 1)
    cell_size_type ncell_max =
        static_cast<cell_size_type>(box->length() / (r_cut_max + r_skin_));
    // determine optimal value from 2) together with b,c)
    size_t warp_size = cuda::device::properties(cuda::device::get()).warp_size();
    double nwarps = particle->nbox / (nu_cell_ * warp_size);
    ncell_ = static_cast<cell_size_type>(ceil(pow(nwarps, 1./dimension)));
    LOG_DEBUG("desired values for number of cells: " << ncell_);
    LOG_DEBUG("upper bound on number of cells: " << ncell_max);
    // respect upper bound
    ncell_ = min(ncell_, ncell_max);

    // compute derived values
    size_t ncells = accumulate(ncell_.begin(), ncell_.end(), 1, multiplies<size_t>());
    cell_size_ = warp_size * static_cast<size_t>(ceil(nwarps / ncells));
    vector_type cell_length_ =
        element_div(static_cast<vector_type>(box->length()), static_cast<vector_type>(ncell_));
    dim_cell_ = cuda::config(
        dim3(
             accumulate(ncell_.begin(), ncell_.end() - 1, 1, multiplies<size_t>())
           , ncell_.back()
        )
      , cell_size_
    );

    if (*min_element(ncell_.begin(), ncell_.end()) < 3) {
        throw std::logic_error("number of cells per dimension must be at least 3");
    }

    LOG("number of placeholders per cell: " << cell_size_);
    LOG("number of cells per dimension: " << ncell_);
    LOG("cell edge lengths: " << cell_length_);
    LOG("desired average cell occupancy: " << nu_cell_);
    double nu_cell_eff = static_cast<double>(particle->nbox) / dim_cell_.threads();
    LOG("effective average cell occupancy: " << nu_cell_eff);

    try {
        cuda::copy(particle->nbox, get_neighbour_kernel<dimension>().nbox);
        cuda::copy(static_cast<fixed_vector<uint, dimension> >(ncell_), get_neighbour_kernel<dimension>().ncell);
        cuda::copy(cell_length_, get_neighbour_kernel<dimension>().cell_length);
        cuda::copy(static_cast<vector_type>(box->length()), get_neighbour_kernel<dimension>().box_length);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw std::logic_error("failed to copy cell parameters to device symbols");
    }

    // volume of n-dimensional sphere with neighbour list radius
    // volume of unit sphere: V_d = π^(d/2) / Γ(1+d/2), Γ(1) = 1, Γ(1/2) = √π
    float_type unit_sphere[5] = {0, 2, M_PI, 4 * M_PI / 3, M_PI * M_PI / 2 };
    assert(dimension <= 4);
    float_type neighbour_sphere = unit_sphere[dimension] * pow(r_cut_max + r_skin_, dimension);
    // number of placeholders per neighbour list
    particle->neighbour_size =
        static_cast<size_t>(ceil(neighbour_sphere * (box->density() / nu_cell_eff)));
    // at least cell_size (or warp_size?) placeholders
    // FIXME what is a sensible lower bound?
    particle->neighbour_size = max(particle->neighbour_size, (unsigned)cell_size_);
    // number of neighbour lists
    particle->neighbour_stride = particle->dim.threads();
    // allocate neighbour lists
    particle->g_neighbour.resize(particle->neighbour_stride * particle->neighbour_size);

    LOG("neighbour list skin: " << r_skin_);
    LOG("number of placeholders per neighbour list: " << particle->neighbour_size);

    try {
        cuda::copy(rr_cut_skin_.data(), g_rr_cut_skin_);
        cuda::copy(particle->neighbour_size, get_neighbour_kernel<dimension>().neighbour_size);
        cuda::copy(particle->neighbour_stride, get_neighbour_kernel<dimension>().neighbour_stride);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw std::logic_error("failed to copy neighbour list parameters to device symbols");
    }

    try {
        g_cell_.resize(dim_cell_.threads());
        g_cell_offset_.resize(dim_cell_.blocks_per_grid());
        g_cell_index_.reserve(particle->dim.threads());
        g_cell_index_.resize(particle->nbox);
        g_cell_permutation_.reserve(particle->dim.threads());
        g_cell_permutation_.resize(particle->nbox);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw std::logic_error("failed to allocate cell placeholders in global device memory");
    }
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void neighbour<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
}

/**
 * Update neighbour lists
 */
template <int dimension, typename float_type>
void neighbour<dimension, float_type>::update()
{
    // rebuild cell lists
    update_cells();
    // rebuild neighbour lists
    update_neighbours();
    // make snapshot of absolute particle displacements
    cuda::copy(particle->g_r, g_r0_);
}

template <int dimension, typename float_type>
typename neighbour_wrapper<dimension>::displacement_impl_type
neighbour<dimension, float_type>::get_displacement_impl(int threads)
{
    switch (threads) {
      case 512:
        return neighbour_wrapper<dimension>::kernel.displacement_impl[0];
      case 256:
        return neighbour_wrapper<dimension>::kernel.displacement_impl[1];
      case 128:
        return neighbour_wrapper<dimension>::kernel.displacement_impl[2];
      case 64:
        return neighbour_wrapper<dimension>::kernel.displacement_impl[3];
      case 32:
        return neighbour_wrapper<dimension>::kernel.displacement_impl[4];
      default:
        throw std::logic_error("invalid reduction thread count");
    }
}

/**
 * Check if neighbour list update is needed
 */
template <int dimension, typename float_type>
bool neighbour<dimension, float_type>::check()
{
    scoped_timer<timer> timer_(at_key<check_>(runtime_));
    try {
        cuda::configure(dim_reduce.grid, dim_reduce.block, dim_reduce.threads_per_block() * sizeof(float));
        displacement_impl(particle->g_r, g_r0_, g_rr_);
        cuda::copy(g_rr_, h_rr_);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw std::logic_error("failed to reduce squared particle displacements on GPU");
    }
    return *max_element(h_rr_.begin(), h_rr_.end()) > rr_skin_half_;
}

/**
 * Update cell lists
 */
template <int dimension, typename float_type>
void neighbour<dimension, float_type>::update_cells()
{
    scoped_timer<timer> timer_(at_key<update_cells_>(runtime_));

    // compute cell indices for particle positions
    cuda::configure(particle->dim.grid, particle->dim.block);
    get_neighbour_kernel<dimension>().compute_cell(particle->g_r, g_cell_index_);

    // generate permutation
    cuda::configure(particle->dim.grid, particle->dim.block);
    get_neighbour_kernel<dimension>().gen_index(g_cell_permutation_);
    sort_(g_cell_index_, g_cell_permutation_);

    // compute global cell offsets in sorted particle list
    cuda::memset(g_cell_offset_, 0xFF);
    cuda::configure(particle->dim.grid, particle->dim.block);
    get_neighbour_kernel<dimension>().find_cell_offset(g_cell_index_, g_cell_offset_);

    // assign particles to cells
    cuda::vector<int> g_ret(1);
    cuda::host::vector<int> h_ret(1);
    cuda::memset(g_ret, EXIT_SUCCESS);
    cuda::configure(dim_cell_.grid, dim_cell_.block);
    get_neighbour_kernel<dimension>().assign_cells(g_ret, g_cell_index_, g_cell_offset_, g_cell_permutation_, g_cell_);
    cuda::copy(g_ret, h_ret);
    if (h_ret.front() != EXIT_SUCCESS) {
        throw std::runtime_error("overcrowded placeholders in cell lists update");
    }
}

/**
 * Update neighbour lists
 */
template <int dimension, typename float_type>
void neighbour<dimension, float_type>::update_neighbours()
{
    scoped_timer<timer> timer_(at_key<update_neighbours_>(runtime_));

    // mark neighbour list placeholders as virtual particles
    cuda::memset(particle->g_neighbour, 0xFF);
    // build neighbour lists
    cuda::vector<int> g_ret(1);
    cuda::host::vector<int> h_ret(1);
    cuda::memset(g_ret, EXIT_SUCCESS);
    cuda::configure(dim_cell_.grid, dim_cell_.block, cell_size_ * (2 + dimension) * sizeof(int));
    get_neighbour_kernel<dimension>().r.bind(particle->g_r);
    get_neighbour_kernel<dimension>().rr_cut_skin.bind(g_rr_cut_skin_);
    get_neighbour_kernel<dimension>().update_neighbours(g_ret, particle->g_neighbour, g_cell_);
    cuda::thread::synchronize();
    cuda::copy(g_ret, h_ret);
    if (h_ret.front() != EXIT_SUCCESS) {
        throw std::runtime_error("overcrowded placeholders in neighbour lists update");
    }
}

template <typename T>
static void register_lua(char const* class_name)
{
    using namespace luabind;
    lua_registry::get()->push_back
    ((
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("gpu")
                [
                    class_<T, shared_ptr<T> >(class_name)
                        .scope
                        [
                            def("options", &T::options)
                        ]
                ]
            ]
        ]
    ));
}

static __attribute__((constructor)) void register_lua()
{
    register_lua<neighbour<3, float> >("neighbour_3_");
    register_lua<neighbour<2, float> >("neighbour_2_");
}

// explicit instantiation
template class neighbour<3, float>;
template class neighbour<2, float>;

}} // namespace mdsim::gpu

template class module<mdsim::gpu::neighbour<3, float> >;
template class module<mdsim::gpu::neighbour<2, float> >;

} // namespace halmd
