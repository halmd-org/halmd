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

#include <algorithm>
#include <boost/bind.hpp>
#include <cmath>
#include <memory>

#include <halmd/mdsim/gpu/integrators/verlet_nvt_andersen.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type, typename RandomNumberGenerator>
verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::verlet_nvt_andersen(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<force_type> force
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<random_type> random
  , float_type timestep
  , float_type temperature
  , float_type coll_rate
  , std::shared_ptr<logger_type> logger
)
  : particle_(particle)
  , force_(force)
  , box_(box)
  , random_(random)
  , coll_rate_(coll_rate)
  , logger_(logger)
{
    set_timestep(timestep);
    set_temperature(temperature);
    LOG("collision rate with heat bath: " << coll_rate_);
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
    sqrt_temperature_ = std::sqrt(temperature_);
    LOG("temperature of heat bath: " << temperature_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::integrate()
{
    cache_proxy<net_force_array_type const> net_force = force_->net_force();
    cache_proxy<position_array_type> position = particle_->position();
    cache_proxy<velocity_array_type> velocity = particle_->velocity();
    cache_proxy<image_array_type> image = particle_->image();

    scoped_timer_type timer(runtime_.integrate);
    try {
        cuda::configure(particle_->dim.grid, particle_->dim.block);
        wrapper_type::kernel.integrate(
            &*position->begin()
          , &*image->begin()
          , &*velocity->begin()
          , &*net_force->begin()
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
    cache_proxy<net_force_array_type const> net_force = force_->net_force();
    cache_proxy<velocity_array_type> velocity = particle_->velocity();

    scoped_timer_type timer(runtime_.finalize);
    try {
        // use CUDA execution dimensions of 'random' since
        // the kernel makes use of the random number generator
        cuda::configure(random_->rng().dim.grid, random_->rng().dim.block);
        wrapper_type::kernel.finalize(
            &*velocity->begin()
          , &*net_force->begin()
          , timestep_
          , sqrt_temperature_
          , coll_prob_
          , particle_->nparticle()
          , particle_->dim.threads()
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
                  , std::shared_ptr<force_type>
                  , std::shared_ptr<box_type const>
                  , std::shared_ptr<random_type>
                  , float_type
                  , float_type
                  , float_type
                  , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_integrators_verlet_nvt_andersen(lua_State* L)
{
    verlet_nvt_andersen<3, float, random::gpu::rand48>::luaopen(L);
    verlet_nvt_andersen<2, float, random::gpu::rand48>::luaopen(L);
    return 0;
}

// explicit instantiation
template class verlet_nvt_andersen<3, float, random::gpu::rand48>;
template class verlet_nvt_andersen<2, float, random::gpu::rand48>;

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
