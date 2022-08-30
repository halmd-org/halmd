/*
 * Copyright © 2019 Roya Ebrahimi Viand
 * Copyright © 2021 Jaslo Ziska
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

#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/geometries/sphere.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace geometries {

template <int dimension, typename float_type>
sphere<dimension, float_type>::sphere(vector_type const& centre, float_type const& radius)
  : centre_(centre)
  , radius_(radius)
  , radius2_(radius_ * radius_)
{}

template <int dimension, typename float_type>
void sphere<dimension, float_type>::log(std::shared_ptr<halmd::logger> logger_) const
{
    LOG("using sphere geometry");
    LOG("centre: " << centre_);
    LOG("radius: " << radius_);
}

template <int dimension, typename float_type>
void sphere<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static std::string class_name("sphere_" + demangled_name<float_type>());
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("geometries")
            [
                class_<sphere>()
              , def(class_name.c_str(), &std::make_shared<sphere
                  , vector_type const&
                  , float_type
                  >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_geometries_sphere(lua_State* L)
{
    sphere<2, float>::luaopen(L);
    sphere<2, double>::luaopen(L);
    sphere<3, float>::luaopen(L);
    sphere<3, double>::luaopen(L);
    return 0;
}

// explicit instantiation
template class sphere<2, float>;
template class sphere<2, double>;
template class sphere<3, float>;
template class sphere<3, double>;

} // namespace geometries
} // namespace mdsim
} // namespace halmd
