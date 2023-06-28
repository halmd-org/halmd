/*
 * Copyright © 2022-2023 Felix Höfling
 * Copyright © 2021      Jaslo Ziska
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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/geometries/cylinder.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <cmath>
#include <string>

namespace halmd {
namespace mdsim {
namespace geometries {

template <int dimension, typename float_type>
cylinder<dimension, float_type>::cylinder(vector_type const& axis, vector_type const& offset, float_type radius)
  : axis_original_(axis)
  , radius_(radius)
  , axis_(axis)
  , offset_(offset)
  , radius2_(radius * radius)
  , length2_4_(length * length / 4)
{
    // normalise axis
    float_type norm = norm_2(axis_);
    if (norm == 0) {
        throw std::invalid_argument("axis of cylinder geometry must be non-zero");
    }
    axis_ /= norm;
}

template <int dimension, typename float_type>
void cylinder<dimension, float_type>::log(std::shared_ptr<halmd::logger> logger_) const
{
    LOG("using cylinder geometry");
    LOG("radius: " << radius_);
    LOG("length: " << length_);
    LOG("axis vector: " << axis_original_);
    LOG_DEBUG("cylinder axis after normalisation: " << axis_);
    LOG("axis offset: " << offset_);
}

template <int dimension, typename float_type>
float_type cuboid<dimension, float_type>::volume() const
{
    return float_type(M_PI) * radius_ * radius_ * length_;
}

template <int dimension, typename float_type>
void cylinder<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static std::string class_name("cylinder_" + demangled_name<float_type>());
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("geometries")
            [
                class_<cylinder>()
                    .property("volume", &cuboid::volume)

              , def(class_name.c_str(), &std::make_shared<cylinder
                  , vector_type const&
                  , vector_type const&
                  , float_type
                  , float_type
                  >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_geometries_cylinder(lua_State* L)
{
    cylinder<2, float>::luaopen(L);
    cylinder<2, double>::luaopen(L);
    cylinder<3, float>::luaopen(L);
    cylinder<3, double>::luaopen(L);
    return 0;
}

// explicit instantiation
template class cylinder<2, float>;
template class cylinder<2, double>;
template class cylinder<3, float>;
template class cylinder<3, double>;

} // namespace geometries
} // namespace mdsim
} // namespace halmd
