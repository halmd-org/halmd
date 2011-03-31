/*
 * Copyright © 2010  Peter Colberg
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

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <luabind/luabind.hpp>
#include <luabind/operator.hpp>

#include <halmd/config.hpp>
#include <halmd/utility/lua/ublas.hpp>

using namespace boost::numeric::ublas;

HALMD_LUA_API int luaopen_libhalmd_ublas(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("ublas")
        [
            class_<vector<double> >("vector_double")
                .def("data", (unbounded_array<double> const& (vector<double>::*)() const) &vector<double>::data)

          , class_<vector<float> >("vector_float")
                .def("data", (unbounded_array<float> const& (vector<float>::*)() const) &vector<float>::data)

          , class_<matrix<double> >("matrix_double")
                .def("data", (unbounded_array<double> const& (matrix<double>::*)() const) &matrix<double>::data)

          , class_<matrix<float> >("matrix_float")
                .def("data", (unbounded_array<float> const& (matrix<float>::*)() const) &matrix<float>::data)

          , class_<symmetric_matrix<double, lower> >("symmetric_matrix_double_lower")
                .def("data", (unbounded_array<double> const& (symmetric_matrix<double, lower>::*)() const) &symmetric_matrix<double, lower>::data)
                .def(tostring(self))
                .def(tostring(const_self))

          , class_<symmetric_matrix<float, lower> >("symmetric_matrix_float_lower")
                .def("data", (unbounded_array<float> const& (symmetric_matrix<float, lower>::*)() const) &symmetric_matrix<float, lower>::data)
                .def(tostring(self))
                .def(tostring(const_self))

          , class_<symmetric_matrix<double, upper> >("symmetric_matrix_double_upper")
                .def("data", (unbounded_array<double> const& (symmetric_matrix<double, upper>::*)() const) &symmetric_matrix<double, upper>::data)
                .def(tostring(self))
                .def(tostring(const_self))

          , class_<symmetric_matrix<float, upper> >("symmetric_matrix_float_upper")
                .def("data", (unbounded_array<float> const& (symmetric_matrix<float, upper>::*)() const) &symmetric_matrix<float, upper>::data)
                .def(tostring(self))
                .def(tostring(const_self))
        ]
    ];
    return 0;
}
