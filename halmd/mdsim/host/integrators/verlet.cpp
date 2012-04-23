/*
 * Copyright © 2008-2012  Peter Colberg
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

#include <halmd/mdsim/host/integrators/verlet.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

template <int dimension, typename float_type>
verlet<dimension, float_type>::verlet(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type const> box
  , double timestep
  , shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
{
    set_timestep(timestep);
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::set_timestep(double timestep)
{
  timestep_ = timestep;
  timestep_half_ = 0.5 * timestep;

  LOG("integration timestep: " << timestep_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::integrate()
{
    scoped_timer_type timer(runtime_.integrate);

    typename particle_type::position_array_type& position = particle_->position();
    typename particle_type::image_array_type& image = particle_->image();
    typename particle_type::velocity_array_type& velocity = particle_->velocity();
    typename particle_type::force_array_type const& force = particle_->force();
    typename particle_type::mass_array_type const& mass = particle_->mass();

    for (size_t i = 0; i < particle_->nparticle(); ++i) {
        vector_type& v = velocity[i] += force[i] * timestep_half_ / mass[i];
        vector_type& r = position[i] += v * timestep_;
        // enforce periodic boundary conditions
        // TODO: reduction is now to (-L/2, L/2) instead of (0, L) as before
        // check that this is OK
        image[i] += box_->reduce_periodic(r);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::finalize()
{
    scoped_timer_type timer(runtime_.finalize);

    typename particle_type::velocity_array_type& velocity = particle_->velocity();
    typename particle_type::force_array_type const& force = particle_->force();
    typename particle_type::mass_array_type const& mass = particle_->mass();

    for (size_t i = 0; i < particle_->nparticle(); ++i) {
        velocity[i] += force[i] * timestep_half_ / mass[i];
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

template <int dimension, typename float_type>
void verlet<dimension, float_type>::luaopen(lua_State* L)
{
    static string const class_name = demangled_name<verlet>();
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("integrators")
                [
                    class_<verlet>(class_name.c_str())
                        .property("integrate", &wrap_integrate<verlet>)
                        .property("finalize", &wrap_finalize<verlet>)
                        .property("timestep", &verlet::timestep)
                        .def("set_timestep", &verlet::set_timestep)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("integrate", &runtime::integrate)
                                .def_readonly("finalize", &runtime::finalize)
                        ]
                        .def_readonly("runtime", &verlet::runtime_)
                ]
            ]

          , namespace_("integrators")
            [
                def("verlet", &make_shared<verlet
                  , shared_ptr<particle_type>
                  , shared_ptr<box_type const>
                  , double
                  , shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_integrators_verlet(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    verlet<3, double>::luaopen(L);
    verlet<2, double>::luaopen(L);
#else
    verlet<3, float>::luaopen(L);
    verlet<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class verlet<3, double>;
template class verlet<2, double>;
#else
template class verlet<3, float>;
template class verlet<2, float>;
#endif

} // namespace mdsim
} // namespace host
} // namespace integrators
} // namespace halmd
