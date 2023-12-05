/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#include <halmd/observables/gpu/density_mode.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
density_mode<dimension, float_type>::density_mode(
    shared_ptr<particle_type const> particle
  , shared_ptr<particle_group_type> particle_group
  , shared_ptr<wavevector_type const> wavevector
  , shared_ptr<logger> logger
)
    // dependency injection
  : particle_(particle)
  , particle_group_(particle_group)
  , wavevector_(wavevector)
  , logger_(logger)
    // member initialisation
  , nq_(wavevector_->value().size())
  , dim_(particle_->dim())
    // memory allocation
  , g_q_(nq_)
  , g_rho_block_(nq_ * dim_.blocks_per_grid())
  , g_rho_(nq_)
  , h_rho_(nq_)
{
    LOG_INFO(
        "CUDA configuration: " << dim_.blocks_per_grid() << " blocks of "
        << dim_.threads_per_block() << " threads each"
    );
    // copy wavevectors to CUDA device
    try {
        // cast from fixed_vector<double, ...> to fixed_vector<float, ...>
        // and finally to gpu_vector_type (float4 or float2)
        auto const& q = wavevector_->value();
        cuda::memory::host::vector<gpu_vector_type> h_q(nq_);
        for (unsigned int i = 0; i < q.size(); ++i) {
            h_q[i] = static_cast<vector_type>(q[i]);
        }
        cuda::copy(h_q.begin(), h_q.end(), g_q_.begin());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy wavevectors to device");
        throw;
    }
}

/**
 * Acquire sample of all density modes from particle group
 */
template <int dimension, typename float_type>
shared_ptr<typename density_mode<dimension, float_type>::result_type const>
density_mode<dimension, float_type>::acquire()
{
    // check validity of caches
    auto const& group_cache  = particle_group_->ordered();
    auto const& position_cache = particle_->position();

    if (group_cache_ != group_cache || position_cache_ != position_cache) {
        // obtain read access to input caches
        auto const& group = read_cache(group_cache);
        auto const& position = read_cache(position_cache);

        LOG_DEBUG("acquire sample");

        scoped_timer_type timer(runtime_.acquire);

        // allocate new memory which allows modules (e.g.,
        // dynamics::blocking_scheme) to hold a previous copy of the result or
        // to track the update via std::weak_ptr.
        result_ = make_shared<result_type>(nq_);

        // compute density modes
        try {
            cuda::texture<gpu_vector_type> t_wavevector(g_q_);

            wrapper_type::kernel.compute.configure(dim_.grid, dim_.block);
            // compute exp(i q·r) for all wavevector/particle pairs and perform block sums
            wrapper_type::kernel.compute(
                t_wavevector
              , position.data(), &*group.begin(), group.size()
              , g_rho_block_.data(), nq_
            );
            cuda::thread::synchronize();

            // finalise block sums for each wavevector
            wrapper_type::kernel.finalise.configure(
                nq_                        // #blocks: one per wavevector, at most 2^31 - 1 ≈ 2 × 10^9
              , dim_.block                 // #threads per block, must be a power of 2
            );
            wrapper_type::kernel.finalise(g_rho_block_.data(), g_rho_.data(), nq_, dim_.blocks_per_grid());
        }
        catch (cuda::error const&) {
            LOG_ERROR("failed to compute density modes on GPU");
            throw;
        }

        // copy data from device and store in density_mode sample
        cuda::copy(g_rho_.begin(), g_rho_.end(), h_rho_.begin());
        auto rho = begin(*result_);
        for (unsigned int i = 0; i < nq_; ++i) {
            // convert from float2 to result type
            *rho++ = fixed_vector<double, 2>({{ h_rho_[i].x, -h_rho_[i].y }}); // FIXME check minus sign on imaginary part
        }

        // update cache observers
        group_cache_ = group_cache;
        position_cache_ = position_cache;
    }

    return result_;
}

template <int dimension, typename float_type>
void density_mode<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                class_<density_mode>()
                    .property("acquisitor", &density_mode::acquisitor)
                    .property("wavevector", &density_mode::wavevector)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("acquire", &runtime::acquire)
                    ]
                    .def_readonly("runtime", &density_mode::runtime_)
            ]
          , def("density_mode", &make_shared<density_mode
              , shared_ptr<particle_type const>
              , shared_ptr<particle_group_type>
              , shared_ptr<wavevector_type const>
              , shared_ptr<logger>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_density_mode(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    density_mode<3, float>::luaopen(L);
    density_mode<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    density_mode<3, dsfloat>::luaopen(L);
    density_mode<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class density_mode<3, float>;
template class density_mode<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class density_mode<3, dsfloat>;
template class density_mode<2, dsfloat>;
#endif

}  // namespace gpu
}  // namespace observables
}  // namespace halmd
