/*
 * Copyright © 2011  Peter Colberg
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

#include <halmd/mdsim/host/positions/excluded_volume.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace positions {

template <int dimension, typename float_type>
excluded_volume<dimension, float_type>::excluded_volume(
    boost::shared_ptr<box_type const> box
  , float_type cell_length
  , boost::shared_ptr<logger_type> logger
)
  : box_(box)
  , logger_(logger)
{
}

template <int dimension, typename float_type>
void excluded_volume<dimension, float_type>::exclude_sphere(
    vector_type const& centre
  , float_type diameter
)
{
}

template <int dimension, typename float_type>
void excluded_volume<dimension, float_type>::exclude_spheres(
    sample_type const& sample
  , std::vector<float_type> diameter
)
{
}

template <int dimension, typename float_type>
bool excluded_volume<dimension, float_type>::place_sphere(
    vector_type const& centre
  , float_type diameter
)
{
}

template <int dimension, typename float_type>
shared_ptr<excluded_volume<dimension, float_type> > make_excluded_volume(
    boost::shared_ptr<typename excluded_volume<dimension, float_type>::box_type const> box
  , float_type cell_length
  , boost::shared_ptr<typename excluded_volume<dimension, float_type>::logger_type> logger
  , boost::shared_ptr<typename excluded_volume<dimension, float_type>::sample_type const>
)
{
    return make_shared<excluded_volume<dimension, float_type> >(box, cell_length, logger);
}

template <int dimension, typename float_type>
void excluded_volume<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name(demangled_name<excluded_volume>());
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("positions")
            [
                class_<excluded_volume>(class_name.c_str())
                    .def("exclude_sphere", &excluded_volume::exclude_sphere)
                    .def("exclude_spheres", &excluded_volume::exclude_spheres)
                    .def("place_sphere", &excluded_volume::place_sphere)

              , def("excluded_volume", &make_excluded_volume<dimension, float_type>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_positions_excluded_volume(lua_State* L)
{
    excluded_volume<3, double>::luaopen(L);
    excluded_volume<2, double>::luaopen(L);
    excluded_volume<3, float>::luaopen(L);
    excluded_volume<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class excluded_volume<3, double>;
template class excluded_volume<2, double>;
template class excluded_volume<3, float>;
template class excluded_volume<2, float>;

} // namespace mdsim
} // namespace host
} // namespace positions
} // namespace halmd
