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

#include <luabind/luabind.hpp>
#include <stdint.h>

#include <halmd/config.hpp>
#include <halmd/utility/signal.hpp>

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
// C++ type name is especially useful for debugging parameter
// mismatches, e.g. if the user tries to register an unsupported
// slot with a signal. Luabind will then print all supported
// slot types, with the exact slot signatures.
//
#define SLOT(...)                                               \
    class_<__VA_ARGS__>(#__VA_ARGS__)                           \
        .def("__call", &__VA_ARGS__::operator())                \

#define CONNECTION(...)                                         \
    class_<__VA_ARGS__>(#__VA_ARGS__)                           \
        .def("disconnect", &__VA_ARGS__::disconnect)            \

/**
 * Lua bindings for halmd::signal<>::slot_function_type.
 *
 * This function registers all slot types used in HALMD for connecting
 * module methods to signals. This allows retrieving a slot from a C++
 * module in Lua, and registering it with a signal in another module,
 * or running the slot directly in Lua.
 */
HALMD_LUA_API int luaopen_libhalmd_utility_lua_signal(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        SLOT( signal<void ()>::slot_function_type )
      , SLOT( signal<void (uint64_t)>::slot_function_type )
      , CONNECTION( signal<void ()>::connection )
      , CONNECTION( signal<void (uint64_t)>::connection )
    ];
    return 0;
}

} // namespace halmd
