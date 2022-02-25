/*
 * Copyright © 2011-2022 Felix Höfling
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

#include <halmd/observables/gpu/kinetic_energy_density_mode.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
kinetic_energy_density_mode<dimension, float_type>::kinetic_energy_density_mode(
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
  , dim_(50, 64 << DEVICE_SCALE) // at most 1024 threads per block
    // memory allocation
  , g_q_(nq_)
  , g_sin_block_(nq_ * dim_.blocks_per_grid()), g_cos_block_(nq_ * dim_.blocks_per_grid())
  , g_sin_(nq_), g_cos_(nq_)
  , h_sin_(nq_), h_cos_(nq_)
{
    LOG_DEBUG(
        "CUDA configuration: " << dim_.blocks_per_grid() << " blocks of "
        << dim_.threads_per_block() << " threads each"
    );
    // copy wavevectors to CUDA device
    try {
        // cast from fixed_vector<double, ...> to fixed_vector<float, ...>
        // and finally to gpu_vector_type (float4 or float2)
        auto const& q = wavevector_->value();
        cuda::host::vector<gpu_vector_type> h_q(nq_);
        for (unsigned int i = 0; i < q.size(); ++i) {
            h_q[i] = static_cast<vector_type>(q[i]);
        }
        cuda::copy(h_q, g_q_);
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
shared_ptr<typename kinetic_energy_density_mode<dimension, float_type>::result_type const>
kinetic_energy_density_mode<dimension, float_type>::acquire()
{
    // check validity of caches
    auto const& group_cache  = particle_group_->ordered();
    auto const& position_cache = particle_->position();
    auto const& velocity_cache = particle_->velocity();

    if (group_cache_ != group_cache || position_cache_ != position_cache || velocity_cache_ != velocity_cache) {
        // obtain read access to input caches
        auto const& group = read_cache(group_cache);
        auto const& position = read_cache(position_cache);
        auto const& velocity = read_cache(velocity_cache);

        LOG_TRACE("acquire sample");

        scoped_timer_type timer(runtime_.acquire);

        // allocate new memory which allows modules (e.g.,
        // dynamics::blocking_scheme) to hold a previous copy of the result or
        // to track the update via std::weak_ptr.
        result_ = make_shared<result_type>(nq_);

        // compute density modes
        try {
            cuda::configure(dim_.grid, dim_.block);
            wrapper_type::kernel.q.bind(g_q_);

            // compute exp(i q·r) for all wavevector/particle pairs and perform block sums
            wrapper_type::kernel.compute(
                position.data(), velocity.data(), &*group.begin(), group.size()
              , g_sin_block_, g_cos_block_, nq_
            );
            cuda::thread::synchronize();

            // finalise block sums for each wavevector
            cuda::configure(
                nq_                        // #blocks: one per wavevector
              , dim_.block                 // #threads per block, must be a power of 2
            );
            wrapper_type::kernel.finalise(g_sin_block_, g_cos_block_, g_sin_, g_cos_, nq_, dim_.blocks_per_grid());
        }
        catch (cuda::error const&) {
            LOG_ERROR("failed to compute density modes on GPU");
            throw;
        }

        // copy data from device and store in density_mode sample
        cuda::copy(g_sin_, h_sin_);
        cuda::copy(g_cos_, h_cos_);
        auto rho_q = begin(*result_);
        for (unsigned int i = 0; i < nq_; ++i) {
            *rho_q++ = {{ h_cos_[i], -h_sin_[i] }};
        }

        // update cache observers
        group_cache_ = group_cache;
        position_cache_ = position_cache;
        velocity_cache_ = velocity_cache;
    }

    return result_;
}

template <int dimension, typename float_type>
void kinetic_energy_density_mode<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                class_<kinetic_energy_density_mode>()
                    .property("acquisitor", &kinetic_energy_density_mode::acquisitor)
                    .property("wavevector", &kinetic_energy_density_mode::wavevector)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("acquire", &runtime::acquire)
                    ]
                    .def_readonly("runtime", &kinetic_energy_density_mode::runtime_)
            ]
          , def("kinetic_energy_density_mode", &make_shared<kinetic_energy_density_mode
              , shared_ptr<particle_type const>
              , shared_ptr<particle_group_type>
              , shared_ptr<wavevector_type const>
              , shared_ptr<logger>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_kinetic_energy_density_mode(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    kinetic_energy_density_mode<3, float>::luaopen(L);
    kinetic_energy_density_mode<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    kinetic_energy_density_mode<3, dsfloat>::luaopen(L);
    kinetic_energy_density_mode<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class kinetic_energy_density_mode<3, float>;
template class kinetic_energy_density_mode<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class kinetic_energy_density_mode<3, dsfloat>;
template class kinetic_energy_density_mode<2, dsfloat>;
#endif

}  // namespace gpu
}  // namespace observables
}  // namespace halmd
