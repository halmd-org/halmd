/*
 * Copyright © 2008-2014  Felix Höfling
 * Copyright © 2013-2014  Nicolas Höft
 * Copyright © 2008-2011  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <boost/lexical_cast.hpp>
#include <exception>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/mdsim/gpu/binning.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

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
 * @param occupancy initial average cell occupancy
 */
template <int dimension, typename float_type>
binning<dimension, float_type>::binning(
    std::shared_ptr<particle_type const> particle
  , std::shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , double skin
  , double occupancy
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
  // allocate parameters
  , r_skin_(skin)
  , r_cut_max_(*std::max_element(r_cut.data().begin(), r_cut.data().end()))
  , device_properties_(cuda::device::get())
{
    LOG("initial cell occupancy: " << occupancy)

    // Find the optimal(?) cell size for given cell occupancy ν.
    // Ideally, we would like to have warp_size placeholders per cell
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
        static_cast<cell_size_type>(box_->length() / (r_cut_max_ + r_skin_));
    // determine optimal value from 2,3) together with b,c)
    size_t warp_size = device_properties_.warp_size();
    double nwarps = particle_->nparticle() / (occupancy * warp_size);
    double volume = std::accumulate(box_->length().begin(), box_->length().end(), 1., std::multiplies<double>());
    ncell_ = static_cast<cell_size_type>(ceil(box_->length() * pow(nwarps / volume, 1./dimension)));
    LOG_DEBUG("desired values for number of cells: " << ncell_);
    LOG_DEBUG("upper bound on number of cells: " << ncell_max);
    // respect upper bound
    ncell_ = element_min(ncell_, ncell_max);

    // compute derived values
    cell_length_ = element_div(static_cast<vector_type>(box_->length()), static_cast<vector_type>(ncell_));

    LOG("neighbour list skin: " << r_skin_);
    LOG("number of cells per dimension: " << ncell_);
    LOG("edge lengths of cells: " << cell_length_);

    size_t ncells = std::accumulate(ncell_.begin(), ncell_.end(), 1, std::multiplies<size_t>());
    set_cell_size(warp_size * static_cast<size_t>(std::ceil(nwarps / ncells)));

    // allocate memory for arrays which are not resized upon changing the
    // number of placeholders
    try {
        g_cell_offset_.resize(dim_cell_.blocks_per_grid());
        g_cell_index_.reserve(particle_->array_size());
        g_cell_index_.resize(particle_->nparticle());
        g_cell_permutation_.reserve(particle_->array_size());
        g_cell_permutation_.resize(particle_->nparticle());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate global device memory");
        throw;
    }
}

template <int dimension, typename float_type>
void binning<dimension, float_type>::set_cell_size(size_t cell_size)
{
    cell_size_ = cell_size;
    size_t ncells = std::accumulate(ncell_.begin(), ncell_.end(), 1, std::multiplies<size_t>());

    LOG("number of placeholders per cell: " << cell_size_);
    LOG("average cell occupancy: " << static_cast<double>(particle_->nparticle()) / (ncells * cell_size_));

    size_t warp_size = device_properties_.warp_size();
    if (cell_size_ % warp_size) {
        throw std::logic_error("cell size must be a multiple of warp size ("
            + boost::lexical_cast<std::string>(warp_size) + ")"
        );
    }

    dim_cell_ = cuda::config(
        dim3(ncells / ncell_.back(), ncell_.back())
      , std::min(cell_size_, device_properties_.max_threads_per_block())
    );
    LOG_DEBUG("CUDA threads per block: " << dim_cell_.threads_per_block());

    try {
        auto g_cell = make_cache_mutable(g_cell_);
        g_cell->resize(ncells * cell_size);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate cell placeholders in global device memory");
        throw;
    }
}

template <int dimension, typename float_type>
cache<typename binning<dimension, float_type>::array_type> const&
binning<dimension, float_type>::g_cell()
{
    auto const& position_cache = particle_->position();
    if (cell_cache_ != position_cache) {
        update();
        cell_cache_ = position_cache;
    }
    return g_cell_;
}

/**
 * Update cell lists
 */
template <int dimension, typename float_type>
void binning<dimension, float_type>::update()
{
    auto const& position = read_cache(particle_->position());
    auto g_cell = make_cache_mutable(g_cell_);

    LOG_TRACE("update cell lists");

    bool overcrowded = false;
    do {
        scoped_timer_type timer(runtime_.update);

        auto const* kernel = &binning_wrapper<dimension>::kernel;
        unsigned int nparticle = particle_->nparticle();

        // compute cell indices for particle positions
        int blockSize = kernel->compute_cell.max_block_size();
        if (!blockSize) blockSize = particle_->dim().block.x;
        cuda::configure(particle_->array_size() / blockSize, blockSize);
        kernel->compute_cell(
            position.data()
          , g_cell_index_
          , cell_length_
          , static_cast<fixed_vector<uint, dimension> >(ncell_)
        );

        // generate permutation
        blockSize = kernel->gen_index.max_block_size();
        if (!blockSize) blockSize = particle_->dim().block.x;
        cuda::configure(particle_->array_size() / blockSize, blockSize);
        kernel->gen_index(g_cell_permutation_, nparticle);
        radix_sort(g_cell_index_.begin(), g_cell_index_.end(), g_cell_permutation_.begin());

        // compute global cell offsets in sorted particle list
        cuda::memset(g_cell_offset_, 0xFF);
        blockSize = kernel->find_cell_offset.max_block_size();
        if (!blockSize) blockSize = particle_->dim().block.x;
        cuda::configure(particle_->array_size() / blockSize, blockSize);
        kernel->find_cell_offset(g_cell_index_, g_cell_offset_, nparticle);

        // assign particles to cells
        cuda::vector<int> g_ret(1);
        cuda::host::vector<int> h_ret(1);
        cuda::memset(g_ret, EXIT_SUCCESS);
        cuda::configure(dim_cell_.grid, dim_cell_.block);
        kernel->assign_cells(
            g_ret
          , g_cell_index_
          , g_cell_offset_
          , g_cell_permutation_
          , &*g_cell->begin()
          , nparticle
          , cell_size_
        );
        cuda::copy(g_ret, h_ret);
        overcrowded = h_ret.front() != EXIT_SUCCESS;
        if (overcrowded) {
            LOG("overcrowded placeholders in cell list update, increase cell size");
            set_cell_size(2 * cell_size_);
        }
    } while(overcrowded);
}

template <int dimension, typename float_type>
void binning<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<binning>()
                .property("r_skin", &binning::r_skin)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("update", &runtime::update)
                ]
                .def_readonly("runtime", &binning::runtime_)

          , def("binning", &std::make_shared<binning
              , std::shared_ptr<particle_type const>
              , std::shared_ptr<box_type const>
              , matrix_type const&
              , double
              , double
              , std::shared_ptr<logger>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_binning(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    binning<3, float>::luaopen(L);
    binning<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    binning<3, dsfloat>::luaopen(L);
    binning<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class binning<3, float>;
template class binning<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class binning<3, dsfloat>;
template class binning<2, dsfloat>;
#endif

} // namespace gpu
} // namespace mdsim
} // namespace halmd
