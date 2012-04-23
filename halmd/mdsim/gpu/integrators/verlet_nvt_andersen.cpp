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
#include <boost/make_shared.hpp>
#include <cmath>

#include <halmd/mdsim/gpu/integrators/verlet_nvt_andersen.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type, typename RandomNumberGenerator>
verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::verlet_nvt_andersen(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type const> box
  , shared_ptr<random_type> random
  , float_type timestep, float_type temperature, float_type coll_rate
  , shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
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
    LOG("integration timestep: " << timestep_);
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::set_temperature(double temperature)
{
    temperature_ = static_cast<float_type>(temperature);
    sqrt_temperature_ = sqrt(temperature_);
    LOG("temperature of heat bath: " << temperature_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::integrate()
{
    try {
        scoped_timer_type timer(runtime_.integrate);
        cuda::configure(
            particle_->dim.grid, particle_->dim.block
        );
        wrapper_type::kernel.integrate(
            particle_->position()
          , particle_->image()
          , particle_->velocity()
          , particle_->force()
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
    // TODO: possibly a performance critical issue:
    // the old implementation had this loop included in update_forces(),
    // which saves one additional read of the forces plus the additional kernel execution
    // and scheduling
    try {
        scoped_timer_type timer(runtime_.finalize);
        // use CUDA execution dimensions of 'random' since
        // the kernel makes use of the random number generator
        cuda::configure(
            random_->rng().dim.grid, random_->rng().dim.block
        );
        wrapper_type::kernel.finalize(
            particle_->velocity()
          , particle_->force()
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

template <typename integrator_type>
static function <void ()>
wrap_integrate(shared_ptr<integrator_type> self)
{
    return bind(&integrator_type::integrate, self);
}

template <typename integrator_type>
static function <void ()>
wrap_finalize(shared_ptr<integrator_type> self)
{
    return bind(&integrator_type::finalize, self);
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::
luaopen(lua_State* L)
{
    static string const class_name = demangled_name<verlet_nvt_andersen>();
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("integrators")
                [
                    class_<verlet_nvt_andersen>(class_name.c_str())
                        .property("integrate", &wrap_integrate<verlet_nvt_andersen>)
                        .property("finalize", &wrap_finalize<verlet_nvt_andersen>)
                        .property("timestep", &verlet_nvt_andersen::timestep)
                        .property("temperature", &verlet_nvt_andersen::temperature)
                        .property("collision_rate", &verlet_nvt_andersen::collision_rate)
                        .def("set_timestep", &verlet_nvt_andersen::set_timestep)
                        .def("set_temperature", &verlet_nvt_andersen::set_temperature)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("integrate", &runtime::integrate)
                                .def_readonly("finalize", &runtime::finalize)
                        ]
                        .def_readonly("runtime", &verlet_nvt_andersen::runtime_)
                ]
            ]

          , namespace_("integrators")
            [
                def("verlet_nvt_andersen", &make_shared<verlet_nvt_andersen
                  , shared_ptr<particle_type>
                  , shared_ptr<box_type const>
                  , shared_ptr<random_type>
                  , float_type
                  , float_type
                  , float_type
                  , shared_ptr<logger_type>
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

} // namespace mdsim
} // namespace gpu
} // namespace integrators
} // namespace halmd
