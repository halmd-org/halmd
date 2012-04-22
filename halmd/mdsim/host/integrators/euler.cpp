/*
 * Copyright © 2011-2012  Michael Kopp and Felix Höfling
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
#include <string>

#include <halmd/mdsim/host/integrators/euler.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

// constructor
template <int dimension, typename float_type>
euler<dimension, float_type>::euler(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type const> box
  , double timestep
  , shared_ptr<logger_type> logger
)
  // dependency injection (initialize public variables)
  : particle_(particle)
  , box_(box)
  , logger_(logger)
{
    set_timestep(timestep);
}

/**
 * set integration timestep
 */
template <int dimension, typename float_type>
void euler<dimension, float_type>::set_timestep(double timestep)
{
  timestep_ = timestep;

  LOG("integration timestep: " << timestep_);
}

/**
 * perform Euler integration: update positions from velocities
 *
 * @f$ r(t + \Delta t) = r(t) + v(t) \Delta t @f$
 */
template <int dimension, typename float_type>
void euler<dimension, float_type>::integrate()
{
    scoped_timer_type timer(runtime_.integrate);

    typename particle_type::position_array_type& position = particle_->position();
    typename particle_type::image_array_type& image = particle_->image();
    typename particle_type::velocity_array_type const& velocity = particle_->velocity();

    for (size_t i = 0 ; i < particle_->nparticle(); ++i)
    {
        vector_type& r = position[i] += velocity[i] * timestep_;
        // enforce periodic boundary conditions
        image[i] += box_->reduce_periodic(r);
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
void euler<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name = demangled_name<euler>();
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("integrators")
                [
                    class_<euler>(class_name.c_str())
                        .property("integrate", &wrap_integrate<euler>)
                        .property("timestep", &euler::set_timestep)
                        .def("set_timestep", &euler::set_timestep)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("integrate", &runtime::integrate)
                        ]
                        .def_readonly("runtime", &euler::runtime_)
                ]
            ]

          , namespace_("integrators")
            [
                def("euler", &make_shared<euler
                  , shared_ptr<particle_type>
                  , shared_ptr<box_type const>
                  , double
                  , shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_integrators_euler(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    euler<3, double>::luaopen(L);
    euler<2, double>::luaopen(L);
#else
    euler<3, float>::luaopen(L);
    euler<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class euler<3, double>;
template class euler<2, double>;
#else
template class euler<3, float>;
template class euler<2, float>;
#endif

} // namespace integrators
} // namespace host
} // namespace mdsim
} // namespace halmd
