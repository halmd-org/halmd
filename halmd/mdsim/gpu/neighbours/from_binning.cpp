/*
 * Copyright © 2008-2015 Felix Höfling
 * Copyright © 2021      Jaslo Ziska
 * Copyright © 2013-2015 Nicolas Höft
 * Copyright © 2008-2011 Peter Colberg
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

#include <halmd/mdsim/gpu/neighbours/from_binning.hpp>
#include <halmd/mdsim/gpu/neighbours/from_binning_kernel.hpp>
#include <halmd/utility/gpu/configure_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/gpu/device.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {

/**
 * construct neighbour list module
 *
 * @param particle mdsim::gpu::particle instances
 * @param binning mdsim::gpu::binning instances
 * @param displacement mdsim::gpu::displacement instances
 * @param box mdsim::box instance
 * @param r_cut force cutoff radius
 * @param skin neighbour list skin
 * @param cell_occupancy desired average cell occupancy
 */
template <int dimension, typename float_type>
from_binning<dimension, float_type>::from_binning(
    std::shared_ptr<particle_type const> particle1
  , std::shared_ptr<particle_type const> particle2
  , std::pair<std::shared_ptr<binning_type>, std::shared_ptr<binning_type>> binning
  , std::pair<std::shared_ptr<displacement_type>, std::shared_ptr<displacement_type>> displacement
  , std::shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , double skin
  , double cell_occupancy
  , std::pair<algorithm, bool> options
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle1_(particle1)
  , particle2_(particle2)
  , binning1_(binning.first)
  , binning2_(binning.second)
  , displacement1_(displacement.first)
  , displacement2_(displacement.second)
  , box_(box)
  , logger_(logger)
  // allocate parameters
  , r_skin_(skin)
  , r_cut_max_(*std::max_element(r_cut.data().begin(), r_cut.data().end()))
  , rr_cut_skin_(r_cut.size1(), r_cut.size2())
  , g_rr_cut_skin_(rr_cut_skin_.data().size())
  , nu_cell_(cell_occupancy) // FIXME neighbour list occupancy
  , preferred_algorithm_(options.first)
  , unroll_force_loop_(options.second)
  , device_properties_(device::get())
{
    for (size_t i = 0; i < r_cut.size1(); ++i) {
        for (size_t j = 0; j < r_cut.size2(); ++j) {
            rr_cut_skin_(i, j) = std::pow(r_cut(i, j) + r_skin_, 2);
        }
    }
    try {
        cuda::copy(rr_cut_skin_.data().begin(), rr_cut_skin_.data().end(),
                   g_rr_cut_skin_.begin());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy neighbour list parameters to device symbols");
        throw;
    }
    set_occupancy(cell_occupancy);
}

template <int dimension, typename float_type>
void from_binning<dimension, float_type>::set_occupancy(double cell_occupancy)
{
    nu_cell_ = cell_occupancy;
    // volume of n-dimensional sphere with neighbour list radius
    // volume of unit sphere: V_d = π^(d/2) / Γ(1+d/2), Γ(1) = 1, Γ(1/2) = √π
    float unit_sphere[5] = {0, 2, M_PI, 4 * M_PI / 3, M_PI * M_PI / 2 };
    assert(dimension <= 4);
    float neighbour_sphere = unit_sphere[dimension] * std::pow(r_cut_max_ + r_skin_, dimension);
    // partial number density
    float density = particle2_->nparticle() / box_->volume();
    // number of placeholders per neighbour list
    size_ = static_cast<size_t>(ceil(neighbour_sphere * (density / cell_occupancy)));
    // at least cell_size (or warp_size?) placeholders
    // FIXME what is a sensible lower bound?
    size_ = std::max(size_, binning2_->cell_size());
    // number of neighbour lists
    stride_ = particle1_->dim().threads();
    // allocate neighbour lists
    auto g_neighbour = make_cache_mutable(g_neighbour_);
    g_neighbour->resize(stride_ * size_);

    LOG("neighbour list skin: " << r_skin_);
    LOG("number of placeholders per neighbour list: " << size_);
}

template <int dimension, typename float_type>
cache<typename from_binning<dimension, float_type>::array_type> const&
from_binning<dimension, float_type>::g_neighbour()
{
    cache<reverse_id_array_type> const& reverse_id_cache1 = particle1_->reverse_id();
    cache<reverse_id_array_type> const& reverse_id_cache2 = particle2_->reverse_id();

    auto current_cache = std::tie(reverse_id_cache1, reverse_id_cache2);

    if (neighbour_cache_ != current_cache || float(displacement1_->compute()) > float(r_skin_ / 2)
        || float(displacement2_->compute()) > float(r_skin_ / 2)) {
        on_prepend_update_();
        update();
        displacement1_->zero();
        displacement2_->zero();
        neighbour_cache_ = current_cache;
        on_append_update_();
    }
    return g_neighbour_;
}

/**
 * Test compatibility of binning parameters with this neighbour list algorithm
 *
 * On the GPU, the binning module is required to have at least 3 cells
 * in each spatial direction in order to be used with the neighbour module.
 */
template <int dimension, typename float_type>
bool from_binning<dimension, float_type>::is_binning_compatible(
    std::shared_ptr<binning_type const> binning1
  , std::shared_ptr<binning_type const> binning2
)
{
    auto ncell = binning2->ncell();
    return *std::min_element(ncell.begin(), ncell.end()) >= 3;
}

/**
 * Update neighbour lists
 */
template <int dimension, typename float_type>
void from_binning<dimension, float_type>::update()
{
    position_array_type const& position1 = read_cache(particle1_->position());
    position_array_type const& position2 = read_cache(particle2_->position());
    cell_array_type const& g_cell2 = read_cache(binning2_->g_cell());
    auto g_neighbour = make_cache_mutable(g_neighbour_);

    LOG_DEBUG("update neighbour lists");

    if (!is_binning_compatible(binning1_, binning2_)) {
        throw std::logic_error("number of cells per dimension must be at least 3");
    }

    // if the number of cells in each spatial direction do not match or
    // the cell sizes are different, the naive implementation is required
    bool use_naive = preferred_algorithm_ == naive
                  || binning1_->ncell() != binning2_->ncell()
                  || binning1_->cell_size() != binning2_->cell_size();

    if (use_naive && preferred_algorithm_ == shared_mem) {
        LOG_WARNING_ONCE("falling back to 'naive' neighbour list algorithm due to incompatible binning modules");
    }

    bool overcrowded = false;
    do {
        scoped_timer_type timer(runtime_.update);

        // mark neighbour list placeholders as virtual particles
        cuda::memset(g_neighbour->begin(), g_neighbour->end(), 0xFF);
        // build neighbour lists
        cuda::memory::device::vector<int> g_ret(1);
        cuda::memory::host::vector<int> h_ret(1);
        cuda::memset(g_ret.begin(), g_ret.end(), EXIT_SUCCESS);

        auto* kernel = &from_binning_wrapper<dimension>::kernel.normal;
        if (unroll_force_loop_) {
            kernel = &from_binning_wrapper<dimension>::kernel.unroll_force_loop;
        }

        if (!use_naive) {
            // update the cell list of binning1 here, as the naive implementation
            // does not need it and calling binning1_->g_cell() may trigger
            // a cell list update
            cell_array_type const& g_cell1 = read_cache(binning1_->g_cell());
            // g_cell() triggers an update and now the cell list sizes may mismatch
            // If so, redo the neighbour list update with the naive implementation
            if (binning1_->cell_size() != binning2_->cell_size()) {
                LOG_WARNING_ONCE("falling back to 'naive' neighbour list algorithm due to mismatching cell list sizes");
                use_naive = true;
                overcrowded = true;    // make sure to retry the neighbour list update
                continue;
            }

            // determine and check size of shared memory
            size_t smem_size = binning2_->cell_size() * (2 + dimension) * sizeof(int);
            if (smem_size > device_properties_.shared_mem_per_block()) {
                LOG_WARNING_ONCE("falling back to 'naive' neighbour list algorithm due to insufficient shared memory");
                use_naive = true;
                overcrowded = true;    // make sure to retry the neighbour list update
                continue;
            }

            cuda::texture<float> rr_cut_skin(g_rr_cut_skin_);
            cuda::texture<float4> r1(position1);
            cuda::texture<float4> r2(position2);

            kernel->update_neighbours.configure(binning2_->dim_cell().grid,
                binning2_->dim_cell().block, smem_size);

            kernel->update_neighbours(
                rr_cut_skin
              , r1
              , r2
              , g_ret
              , &*g_neighbour->begin()
              , size_
              , stride_
              , &*g_cell1.begin()
              , &*g_cell2.begin()
              , rr_cut_skin_.size1()
              , rr_cut_skin_.size2()
              , g_cell1.size()
              , binning2_->ncell()
              , static_cast<vector_type>(box_->length())
            );
        }
        else {
            configure_kernel(kernel->update_neighbours_naive, particle1_->dim(),
                false);

            cuda::texture<float> rr_cut_skin(g_rr_cut_skin_);
            cuda::texture<float4> r2(position2);

            kernel->update_neighbours_naive(
                rr_cut_skin
              , r2
              , g_ret
              , position1.data()
              , particle1_->nparticle()
              , particle1_ == particle2_
              , &*g_neighbour->begin()
              , size_
              , stride_
              , &*g_cell2.begin()
              , rr_cut_skin_.size1()
              , rr_cut_skin_.size2()
              , binning2_->ncell()
              , binning2_->cell_length()
              , binning2_->cell_size()
              , static_cast<vector_type>(box_->length())
            );
        }
        cuda::thread::synchronize();
        cuda::copy(g_ret.begin(), g_ret.end(), h_ret.begin());
        overcrowded = h_ret.front() != EXIT_SUCCESS;
        if (overcrowded) {
            LOG("overcrowded placeholders in neighbour lists update, reducing occupancy");
            set_occupancy(nu_cell_ / 2);
        }
    } while (overcrowded);
}

template <int dimension, typename float_type>
float from_binning<dimension, float_type>::defaults::occupancy() {
    return 0.4;
}

template<typename float_type>
struct variant_name;

template<>
struct variant_name<float>
{
    static constexpr const char *name = "float";
};

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template<>
struct variant_name<dsfloat>
{
    static constexpr const char *name = "dsfloat";
};
#endif

template <int dimension, typename float_type>
void from_binning<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    std::string const defaults_name("defaults_" + std::to_string(dimension) + "_" + std::string(variant_name<float_type>::name));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("neighbours")
            [
                class_<from_binning, _Base>()
                    .property("r_skin", &from_binning::r_skin)
                    .property("cell_occupancy", &from_binning::cell_occupancy)
                    .def("on_prepend_update", &from_binning::on_prepend_update)
                    .def("on_append_update", &from_binning::on_append_update)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("update", &runtime::update)
                    ]
                    .def_readonly("runtime", &from_binning::runtime_)
              , def("from_binning", &std::make_shared<from_binning
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<particle_type const>
                  , std::pair<std::shared_ptr<binning_type>, std::shared_ptr<binning_type>>
                  , std::pair<std::shared_ptr<displacement_type>, std::shared_ptr<displacement_type>>
                  , std::shared_ptr<box_type const>
                  , matrix_type const&
                  , double
                  , double
                  , std::pair<algorithm, bool>
                  , std::shared_ptr<logger>
                >)
              , def("is_binning_compatible", &from_binning::is_binning_compatible)
            ]
          , namespace_(defaults_name.c_str())
            [
                namespace_("from_binning")
                [
                    def("occupancy", &defaults::occupancy)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_neighbours_from_binning(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    from_binning<3, float>::luaopen(L);
    from_binning<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    from_binning<3, dsfloat>::luaopen(L);
    from_binning<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class from_binning<3, float>;
template class from_binning<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class from_binning<3, dsfloat>;
template class from_binning<2, dsfloat>;
#endif

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd
