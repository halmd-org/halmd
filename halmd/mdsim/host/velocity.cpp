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

#include <boost/foreach.hpp>

#include <halmd/mdsim/host/velocity.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
velocity<dimension, float_type>::velocity(
    boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , logger_(logger)
{
}

/**
 * Rescale magnitude of all velocities by 'factor'
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::rescale(double factor)
{
    BOOST_FOREACH (vector_type& v, particle_->v) {
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
    BOOST_FOREACH (vector_type& v, particle_->v) {
        v += delta;
    }
}

/**
 * First shift, then rescale all velocities
 */
template <int dimension, typename float_type>
void velocity<dimension, float_type>::shift_rescale(vector_type const& delta, double factor)
{
    BOOST_FOREACH (vector_type& v, particle_->v) {
        v += delta;
        v *= factor;
    }
}

template <int dimension, typename float_type>
void velocity<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("velocity_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                class_<velocity, boost::shared_ptr<_Base>, bases<_Base> >(class_name.c_str())
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_velocity(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    velocity<3, double>::luaopen(L);
    velocity<2, double>::luaopen(L);
#else
    velocity<3, float>::luaopen(L);
    velocity<2, float>::luaopen(L);
#endif
    return 0;
}

#ifndef USE_HOST_SINGLE_PRECISION
template class velocity<3, double>;
template class velocity<2, double>;
#else
template class velocity<3, float>;
template class velocity<2, float>;
#endif

} // namespace mdsim
} // namespace host
} // namespace halmd
