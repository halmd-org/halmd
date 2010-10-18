/*
 * Copyright © 2010  Peter Colberg and Felix Höfling
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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/velocity.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host
{

template <int dimension, typename float_type>
velocity<dimension, float_type>::velocity(
    shared_ptr<particle_type> particle
)
  // dependency injection
  : particle(particle)
{
}

/**
 * Rescale magnitude of all velocities by 'factor'
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::rescale(double factor)
{
    BOOST_FOREACH (vector_type& v, particle->v) {
        v *= factor;
    }
    LOG("velocities rescaled by factor " << factor);
}

/**
 * Shift all velocities by 'delta'
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::shift(vector_type const& delta)
{
    BOOST_FOREACH (vector_type& v, particle->v) {
        v += delta;
    }
}

/**
 * First shift, then rescale all velocities
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::shift_rescale(vector_type const& delta, double factor)
{
    BOOST_FOREACH (vector_type& v, particle->v) {
        v += delta;
        v *= factor;
    }
}

template <int dimension, typename float_type>
void velocity<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    string class_name("velocity_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    class_<velocity, shared_ptr<_Base>, bases<_Base> >(class_name.c_str())
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &velocity<3, double>::luaopen
    ]
    [
        &velocity<2, double>::luaopen
    ];
#else
    [
        &velocity<3, float>::luaopen
    ]
    [
        &velocity<2, float>::luaopen
    ];
#endif
}

#ifndef USE_HOST_SINGLE_PRECISION
template class velocity<3, double>;
template class velocity<2, double>;
#else
template class velocity<3, float>;
template class velocity<2, float>;
#endif

}} // namespace mdsim::host

} // namespace halmd
