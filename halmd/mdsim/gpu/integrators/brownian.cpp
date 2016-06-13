/*
 * Copyright © 2011-2012  Michael Kopp and Felix Höfling
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

#include <memory>
#include <iostream>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/integrators/brownian.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type, typename RandomNumberGenerator>
brownian<dimension, float_type, RandomNumberGenerator>::brownian(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<random_type> random
  , std::shared_ptr<box_type const> box
  , double timestep
  , double T
  , matrix_type const& D
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , random_(random)
  , box_(box)
  , temperature_(T)
  , D_(D)
  , g_param_(D_.size1() * D_.size2() )
  , logger_(logger)
{
    set_timestep(timestep);
    
    cuda::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p;
        p[0] = D_(i,0);
        p[1] = D_(i,1);
        p[2] = D_(i,2);
        p[3] = D_(i,3);
        param[i] = p;
    } 

    cuda::copy(param, g_param_);
}

/**
 * set integration timestep
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void brownian<dimension, float_type, RandomNumberGenerator>::set_timestep(double timestep)
{
    timestep_ = timestep;
}

/**
 * set temperature of the heat bath
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void brownian<dimension, float_type, RandomNumberGenerator>::set_temperature(double temperature)
{
    temperature_ = temperature;
}
/**
 * perform Brownian integration: update positions from random distribution 
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void brownian<dimension, float_type, RandomNumberGenerator>::integrate()
{
    LOG_TRACE("update positions")

    velocity_array_type const& velocity = read_cache(particle_->velocity());
    force_array_type const& force = read_cache(particle_->force());
    torque_array_type const& torque = read_cache(particle_->torque());

    // invalidate the particle caches after accessing the velocity!
    auto position = make_cache_mutable(particle_->position());
    auto orientation = make_cache_mutable(particle_->orientation());
    auto image = make_cache_mutable(particle_->image());
    bind_textures();
    scoped_timer_type timer(runtime_.integrate);

    try {
        // use CUDA execution dimensions of 'random' since
        // the kernel makes use of the random number generator
        cuda::configure(random_->rng().dim.grid, random_->rng().dim.block);
        wrapper_type::kernel.integrate(
            &*position->begin()
          , &*orientation->begin()
          , &*image->begin()
          , &*velocity.begin()
          , &*force.begin()
          , &*torque.begin()
          , timestep_
          , temperature_
          , random_->rng().rng()
          , particle_->nparticle()
          , particle_->dim.threads()
          , static_cast<vector_type>(box_->length())
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream Brownian integration on GPU");
        throw;
    }
}

template <typename integrator_type>
static std::function<void ()>
wrap_integrate(std::shared_ptr<integrator_type> self)
{
    return [=]() {
        self->integrate();
    };
}

template <typename integrator_type>
static std::function<void ()>
wrap_finalize(std::shared_ptr<integrator_type> self)
{
    return [=]() {
        self->finalize();
    };
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void brownian<dimension, float_type, RandomNumberGenerator>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<brownian>()
                    .property("integrate", &wrap_integrate<brownian>)
                    .property("timestep", &brownian::timestep)
                    .property("temperature", &brownian::temperature)
                    .def("set_timestep", &brownian::set_timestep)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("integrate", &runtime::integrate)
                    ]
                    .def_readonly("runtime", &brownian::runtime_)

              , def("brownian", &std::make_shared<brownian
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<random_type>
                  , std::shared_ptr<box_type const>
                  , double
                  , double
                  , matrix_type const&
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_integrators_brownian(lua_State* L)
{
    brownian<3, float, halmd::random::gpu::rand48>::luaopen(L);
    brownian<3, float, halmd::random::gpu::mrg32k3a>::luaopen(L);
    //brownian<2, float, halmd::random::gpu::rand48>::luaopen(L);
    return 0;
}

// explicit instantiation
template class brownian<3, float, halmd::random::gpu::rand48>;
template class brownian<3, float, halmd::random::gpu::mrg32k3a>;
//let's stay 3d for now
//template class brownian<2, float, halmd::random::gpu::rand48>;

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
