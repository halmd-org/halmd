/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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
#include <boost/shared_ptr.hpp>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/forces/zero.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace forces {

template <int dimension, typename float_type>
zero<dimension, float_type>::zero(shared_ptr<particle_type> particle)
  // dependency injection
  : particle(particle)
  // memory allocation
  , stress_pot_(particle->nbox)
{
}

template <int dimension, typename float_type>
void zero<dimension, float_type>::compute()
{
    LOG_TRACE("zero particle forces");
    fill(particle->f.begin(), particle->f.end(), 0);
    fill(stress_pot_.begin(), stress_pot_.end(), 0); // FIXME check aux flag
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(zero<dimension, float_type> const&)
{
    return zero<dimension, float_type>::module_name();
}

template <int dimension, typename float_type>
void zero<dimension, float_type>::luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static string class_name("zero_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("forces")
                [
                    class_<zero, shared_ptr<_Base_Base>, bases<_Base_Base, _Base> >(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type>
                        >())
                        .property("module_name", &module_name_wrapper<dimension, float_type>)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_forces_zero(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    zero<3, double>::luaopen(L);
    zero<2, double>::luaopen(L);
#else
    zero<3, float>::luaopen(L);
    zero<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class zero<3, double>;
template class zero<2, double>;
#else
template class zero<3, float>;
template class zero<2, float>;
#endif

} // namespace mdsim
} // namespace host
} // namespace forces
} // namespace halmd
