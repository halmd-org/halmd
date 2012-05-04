/*
 * Copyright © 2011  Felix Höfling
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

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <string>

#include <halmd/observables/density_mode.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {

template <typename density_mode_type>
typename density_mode_type::slot_function_type
acquire_wrapper(boost::shared_ptr<density_mode_type> density_mode)
{
    return bind(&density_mode_type::acquire, density_mode);
}

template <typename density_mode_type>
static function<vector<double> const& ()>
wrap_wavenumber(boost::shared_ptr<density_mode_type const> density_mode)
{
    return bind(&density_mode_type::wavenumber, density_mode);
}

template <int dimension>
void density_mode<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("density_mode_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<density_mode, boost::shared_ptr<density_mode> >(class_name.c_str())
                .property("value", &density_mode::value)
                .property("wavenumber", &wrap_wavenumber<density_mode>)
                .property("acquire", &acquire_wrapper<density_mode>)
                .def("on_acquire", &density_mode::on_acquire)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_density_mode(lua_State* L)
{
    density_mode<3>::luaopen(L);
    density_mode<2>::luaopen(L);
    return 0;
}

} // namespace observables
} // namespace halmd
