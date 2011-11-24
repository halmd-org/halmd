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

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/function.hpp>
#include <luabind/luabind.hpp>
#include <stdint.h> // uint64_t
#include <vector>

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/lua/fixed_vector_converter.hpp>
#include <halmd/utility/lua/vector_converter.hpp>
#include <halmd/utility/raw_allocator.hpp>

using namespace boost;
using namespace std;

namespace halmd {

// This macro uses __VA_ARGS__ to support template types with commas.
// __VA_ARGS__ is part of the C99 standard, and will be part of C++0x,
// therefore most C++ compilers already support __VA_ARGS__ as an
// extension. This was tested with GCC 4.4 and Clang 2.9.
//
// The stringification turns the C++ type name into a Lua class name.
// A Lua class name may be any string of characters, e.g. spaces,
// commas, brackets or ampersands. As the registered classes are
// never constructed in Lua, but returned from C++ modules, the
// class names only have informational purposes. Use of the full
// C++ type name is especially useful for debugging argument
// mismatches, e.g. if the user tries to register an unsupported
// data slot with the H5MD writer. Luabind will then print all
// supported slot types, with the exact slot signatures.
//
#define SLOT(...)                                               \
    class_<__VA_ARGS__>(#__VA_ARGS__)                           \
        .def("__call", &__VA_ARGS__::operator())                \

/*
 * Lua bindings for boost::function<> with return value.
 *
 * This function registers all data slot types used in HALMD.
 * This allows retrieving a data slot from a C++ module in Lua, and
 * registering it with the H5MD or other writer, or running the slot
 * directly in Lua to convert the data to a Lua number or table.
 */
HALMD_LUA_API int luaopen_libhalmd_utility_lua_function(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        SLOT( function<float ()> )
      , SLOT( function<float& ()> )
      , SLOT( function<float const& ()> )
      , SLOT( function<double ()> )
      , SLOT( function<double& ()> )
      , SLOT( function<double const& ()> )
      , SLOT( function<fixed_vector<float, 2> ()> )
      , SLOT( function<fixed_vector<float, 2>& ()> )
      , SLOT( function<fixed_vector<float, 2> const& ()> )
      , SLOT( function<fixed_vector<float, 3> ()> )
      , SLOT( function<fixed_vector<float, 3>& ()> )
      , SLOT( function<fixed_vector<float, 3> const& ()> )
      , SLOT( function<fixed_vector<double, 2> ()> )
      , SLOT( function<fixed_vector<double, 2>& ()> )
      , SLOT( function<fixed_vector<double, 2> const& ()> )
      , SLOT( function<fixed_vector<double, 3> ()> )
      , SLOT( function<fixed_vector<double, 3>& ()> )
      , SLOT( function<fixed_vector<double, 3> const& ()> )
      , SLOT( function<vector<float> ()> )
      , SLOT( function<vector<float>& ()> )
      , SLOT( function<vector<float> const& ()> )
      , SLOT( function<vector<double> ()> )
      , SLOT( function<vector<double>& ()> )
      , SLOT( function<vector<double> const& ()> )
      , SLOT( function<vector<fixed_vector<float, 2> > ()> )
      , SLOT( function<vector<fixed_vector<float, 2> >& ()> )
      , SLOT( function<vector<fixed_vector<float, 2> > const& ()> )
      , SLOT( function<vector<fixed_vector<float, 3> > ()> )
      , SLOT( function<vector<fixed_vector<float, 3> >& ()> )
      , SLOT( function<vector<fixed_vector<float, 3> > const& ()> )
      , SLOT( function<vector<fixed_vector<double, 2> > ()> )
      , SLOT( function<vector<fixed_vector<double, 2> >& ()> )
      , SLOT( function<vector<fixed_vector<double, 2> > const& ()> )
      , SLOT( function<vector<fixed_vector<double, 3> > ()> )
      , SLOT( function<vector<fixed_vector<double, 3> >& ()> )
      , SLOT( function<vector<fixed_vector<double, 3> > const& ()> )
      , SLOT( function<vector<fixed_vector<float, 2>, raw_allocator<fixed_vector<float, 2> > >& ()> )
      , SLOT( function<vector<fixed_vector<float, 2>, raw_allocator<fixed_vector<float, 2> > > const& ()> )
      , SLOT( function<vector<fixed_vector<float, 3>, raw_allocator<fixed_vector<float, 3> > >& ()> )
      , SLOT( function<vector<fixed_vector<float, 3>, raw_allocator<fixed_vector<float, 3> > > const& ()> )
      , SLOT( function<vector<fixed_vector<double, 2>, raw_allocator<fixed_vector<double, 2> > >& ()> )
      , SLOT( function<vector<fixed_vector<double, 2>, raw_allocator<fixed_vector<double, 2> > > const& ()> )
      , SLOT( function<vector<fixed_vector<double, 3>, raw_allocator<fixed_vector<double, 3> > >& ()> )
      , SLOT( function<vector<fixed_vector<double, 3>, raw_allocator<fixed_vector<double, 3> > > const& ()> )
      , SLOT( function<vector<unsigned int, raw_allocator<unsigned int > >& ()> )
      , SLOT( function<vector<unsigned int, raw_allocator<unsigned int > > const& ()> )
      , SLOT( function<vector<array<float, 3> > ()> )
      , SLOT( function<vector<array<float, 3> >& ()> )
      , SLOT( function<vector<array<float, 3> > const& ()> )
      , SLOT( function<vector<array<double, 3> > ()> )
      , SLOT( function<vector<array<double, 3> >& ()> )
      , SLOT( function<vector<array<double, 3> > const& ()> )
      , SLOT( function<multi_array<float, 2> ()> )
      , SLOT( function<multi_array<float, 2>& ()> )
      , SLOT( function<multi_array<float, 2> const& ()> )
      , SLOT( function<multi_array<float, 3> ()> )
      , SLOT( function<multi_array<float, 3>& ()> )
      , SLOT( function<multi_array<float, 3> const& ()> )
      , SLOT( function<multi_array<double, 2> ()> )
      , SLOT( function<multi_array<double, 2>& ()> )
      , SLOT( function<multi_array<double, 2> const& ()> )
      , SLOT( function<multi_array<double, 3> ()> )
      , SLOT( function<multi_array<double, 3>& ()> )
      , SLOT( function<multi_array<double, 3> const& ()> )
      , SLOT( function<multi_array<uint64_t, 2> ()> )
      , SLOT( function<multi_array<uint64_t, 2>& ()> )
      , SLOT( function<multi_array<uint64_t, 2> const& ()> )
      , SLOT( function<multi_array<uint64_t, 3> ()> )
      , SLOT( function<multi_array<uint64_t, 3>& ()> )
      , SLOT( function<multi_array<uint64_t, 3> const& ()> )
    ];
    return 0;
}

} // namespace halmd
