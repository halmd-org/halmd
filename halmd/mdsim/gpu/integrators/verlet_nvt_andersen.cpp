/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <algorithm>
#include <boost/bind/bind.hpp>
#include <cmath>
#include <memory>

#include <halmd/mdsim/gpu/integrators/verlet_nvt_andersen.hpp>
#include <halmd/utility/gpu/configure_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type, typename RandomNumberGenerator>
verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::verlet_nvt_andersen(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<random_type> random
  , double timestep
  , double temperature
  , double coll_rate
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , box_(box)
  , random_(random)
  , coll_rate_(coll_rate)
  , logger_(logger)
{
    set_timestep(timestep);
    set_temperature(temperature);
    LOG("collision rate with heat bath: " << float(coll_rate_));
}

/**
 * set integration time-step
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::set_timestep(double timestep)
{
    timestep_ = timestep;
    coll_prob_ = coll_rate_ * timestep;
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::set_temperature(double temperature)
{
    temperature_ = temperature;
    sqrt_temperature_ = std::sqrt(temperature);
    LOG("temperature of heat bath: " << float(temperature_));
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::integrate()
{
    force_array_type const& force = read_cache(particle_->force());

    LOG_TRACE("update positions and velocities: first leapfrog half-step");
    scoped_timer_type timer(runtime_.integrate);

    // invalidate the particle caches after accessing the force!
    auto position = make_cache_mutable(particle_->position());
    auto velocity = make_cache_mutable(particle_->velocity());
    auto image = make_cache_mutable(particle_->image());

    try {
        configure_kernel(wrapper_type::kernel.integrate, particle_->dim(), true);
        wrapper_type::kernel.integrate(
            position->data()
          , image->data()
          , velocity->data()
          , force.data()
          , timestep_
          , static_cast<vector_type>(box_->length())
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream first leapfrog step on GPU");
        throw;
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::finalize()
{
    force_array_type const& force = read_cache(particle_->force());

    LOG_TRACE("update velocities: second leapfrog half-step");
    scoped_timer_type timer(runtime_.finalize);

    // invalidate the particle caches after accessing the force!
    auto velocity = make_cache_mutable(particle_->velocity());

    try {
        // use CUDA execution dimensions of 'random' since
        // the kernel makes use of the random number generator
        wrapper_type::kernel.finalize.configure(random_->rng().dim.grid,
            random_->rng().dim.block);
        wrapper_type::kernel.finalize(
            velocity->data()
          , force.data()
          , timestep_
          , sqrt_temperature_
          , coll_prob_
          , particle_->nparticle()
          , particle_->dim().threads()
          , random_->rng().rng()
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream second leapfrog step on GPU");
        throw;
    }
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<verlet_nvt_andersen>()
                    .def("integrate", &verlet_nvt_andersen::integrate)
                    .def("finalize", &verlet_nvt_andersen::finalize)
                    .def("set_timestep", &verlet_nvt_andersen::set_timestep)
                    .def("set_temperature", &verlet_nvt_andersen::set_temperature)
                    .property("timestep", &verlet_nvt_andersen::timestep)
                    .property("temperature", &verlet_nvt_andersen::temperature)
                    .property("collision_rate", &verlet_nvt_andersen::collision_rate)
                    .scope
                    [
                        class_<runtime>()
                            .def_readonly("integrate", &runtime::integrate)
                            .def_readonly("finalize", &runtime::finalize)
                    ]
                    .def_readonly("runtime", &verlet_nvt_andersen::runtime_)

              , def("verlet_nvt_andersen", &std::make_shared<verlet_nvt_andersen
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type const>
                  , std::shared_ptr<random_type>
                  , double
                  , double
                  , double
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_integrators_verlet_nvt_andersen(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    verlet_nvt_andersen<3, float, random::gpu::mrg32k3a>::luaopen(L);
    verlet_nvt_andersen<2, float, random::gpu::mrg32k3a>::luaopen(L);
    verlet_nvt_andersen<3, float, random::gpu::rand48>::luaopen(L);
    verlet_nvt_andersen<2, float, random::gpu::rand48>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    verlet_nvt_andersen<3, dsfloat, random::gpu::mrg32k3a>::luaopen(L);
    verlet_nvt_andersen<2, dsfloat, random::gpu::mrg32k3a>::luaopen(L);
    verlet_nvt_andersen<3, dsfloat, random::gpu::rand48>::luaopen(L);
    verlet_nvt_andersen<2, dsfloat, random::gpu::rand48>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class verlet_nvt_andersen<3, float, random::gpu::mrg32k3a>;
template class verlet_nvt_andersen<2, float, random::gpu::mrg32k3a>;
template class verlet_nvt_andersen<3, float, random::gpu::rand48>;
template class verlet_nvt_andersen<2, float, random::gpu::rand48>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class verlet_nvt_andersen<3, dsfloat, random::gpu::mrg32k3a>;
template class verlet_nvt_andersen<2, dsfloat, random::gpu::mrg32k3a>;
template class verlet_nvt_andersen<3, dsfloat, random::gpu::rand48>;
template class verlet_nvt_andersen<2, dsfloat, random::gpu::rand48>;
#endif

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
