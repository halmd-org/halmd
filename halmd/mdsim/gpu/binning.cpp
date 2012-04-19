/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <boost/bind.hpp>
#include <exception>

#include <halmd/mdsim/gpu/binning.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * construct particle binning module
 *
 * @param particle mdsim::gpu::particle instance
 * @param box mdsim::box instance
 * @param cutoff force cutoff radius
 * @param skin neighbour list skin
 * @param cell_occupancy desired average cell occupancy
 */
template <int dimension, typename float_type>
binning<dimension, float_type>::binning(
    shared_ptr<particle_type const> particle
  , shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , double skin
  , shared_ptr<logger_type> logger
  , double cell_occupancy
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
  // allocate parameters
  , r_skin_(skin)
  , nu_cell_(cell_occupancy)
  , sort_(particle_->nparticle(), particle_->dim.threads_per_block())
{
    float_type r_cut_max = *max_element(r_cut.data().begin(), r_cut.data().end());
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
    // 3) n_cell respects the aspect ratios of the simulation box
    // 4) ν_eff ≤ ν  (actual vs. desired occupancy)
    //
    // upper bound on ncell_ provided by 1)
    cell_size_type ncell_max =
        static_cast<cell_size_type>(box_->length() / (r_cut_max + r_skin_));
    // determine optimal value from 2,3) together with b,c)
    size_t warp_size = cuda::device::properties(cuda::device::get()).warp_size();
    double nwarps = particle_->nparticle() / (nu_cell_ * warp_size);
    double volume = accumulate(box_->length().begin(), box_->length().end(), 1., multiplies<double>());
    ncell_ = static_cast<cell_size_type>(ceil(box_->length() * pow(nwarps / volume, 1./dimension)));
    LOG_DEBUG("desired values for number of cells: " << ncell_);
    LOG_DEBUG("upper bound on number of cells: " << ncell_max);
    // respect upper bound
    ncell_ = element_min(ncell_, ncell_max);

    // compute derived values
    size_t ncells = accumulate(ncell_.begin(), ncell_.end(), 1, multiplies<size_t>());
    cell_size_ = warp_size * static_cast<size_t>(ceil(nwarps / ncells));
    cell_length_ = element_div(static_cast<vector_type>(box_->length()), static_cast<vector_type>(ncell_));
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
    nu_cell_eff_ = static_cast<double>(particle_->nparticle()) / dim_cell_.threads();
    LOG("effective average cell occupancy: " << nu_cell_eff_);

    try {
        cuda::copy(particle_->nparticle(), get_binning_kernel<dimension>().nbox);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy cell parameters to device symbols");
        throw;
    }

    try {
        g_cell_.resize(dim_cell_.threads());
        g_cell_offset_.resize(dim_cell_.blocks_per_grid());
        g_cell_index_.reserve(particle_->dim.threads());
        g_cell_index_.resize(particle_->nparticle());
        g_cell_permutation_.reserve(particle_->dim.threads());
        g_cell_permutation_.resize(particle_->nparticle());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate cell placeholders in global device memory");
        throw;
    }
}

/**
 * Update cell lists
 */
template <int dimension, typename float_type>
void binning<dimension, float_type>::update()
{
    LOG_TRACE("update cell lists");

    scoped_timer_type timer(runtime_.update);

    // compute cell indices for particle positions
    cuda::configure(particle_->dim.grid, particle_->dim.block);
    get_binning_kernel<dimension>().compute_cell(
        particle_->g_r
      , g_cell_index_
      , cell_length_
      , static_cast<fixed_vector<uint, dimension> >(ncell_)
    );

    // generate permutation
    cuda::configure(particle_->dim.grid, particle_->dim.block);
    get_binning_kernel<dimension>().gen_index(g_cell_permutation_);
    sort_(g_cell_index_, g_cell_permutation_);

    // compute global cell offsets in sorted particle list
    cuda::memset(g_cell_offset_, 0xFF);
    cuda::configure(particle_->dim.grid, particle_->dim.block);
    get_binning_kernel<dimension>().find_cell_offset(g_cell_index_, g_cell_offset_);

    // assign particles to cells
    cuda::vector<int> g_ret(1);
    cuda::host::vector<int> h_ret(1);
    cuda::memset(g_ret, EXIT_SUCCESS);
    cuda::configure(dim_cell_.grid, dim_cell_.block);
    get_binning_kernel<dimension>().assign_cells(g_ret, g_cell_index_, g_cell_offset_, g_cell_permutation_, g_cell_);
    cuda::copy(g_ret, h_ret);
    if (h_ret.front() != EXIT_SUCCESS) {
        throw std::runtime_error("overcrowded placeholders in cell lists update");
    }
}

template <typename binning_type>
typename signal<void ()>::slot_function_type
wrap_update(shared_ptr<binning_type> binning)
{
    return bind(&binning_type::update, binning);
}

template <int dimension, typename float_type>
float_type binning<dimension, float_type>::defaults::occupancy() {
    return 0.4;
}

template <int dimension, typename float_type>
void binning<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("binning_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                class_<binning, shared_ptr<binning> >(class_name.c_str())
                    .def(constructor<
                        shared_ptr<particle_type const>
                      , shared_ptr<box_type const>
                      , matrix_type const&
                      , double
                      , shared_ptr<logger_type>
                      , double
                    >())
                    .property("update", &wrap_update<binning>)
                    .property("r_skin", &binning::r_skin)
                    .property("cell_occupancy", &binning::cell_occupancy)
                    .scope
                    [
                        namespace_("defaults")
                        [
                            def("occupancy", &defaults::occupancy)
                        ]

                      , class_<runtime>("runtime")
                            .def_readonly("update", &runtime::update)
                    ]
                    .def_readonly("runtime", &binning::runtime_)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_binning(lua_State* L)
{
    binning<3, float>::luaopen(L);
    binning<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class binning<3, float>;
template class binning<2, float>;

} // namespace mdsim
} // namespace gpu
} // namespace halmd
