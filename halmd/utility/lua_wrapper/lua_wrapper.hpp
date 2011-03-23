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

#ifndef HALMD_UTILITY_LUA_WRAPPER_LUA_WRAPPER_HPP
#define HALMD_UTILITY_LUA_WRAPPER_LUA_WRAPPER_HPP

#include <boost/shared_ptr.hpp> // pointer holder for luabind::class_
#include <lua.hpp>
#include <luabind/luabind.hpp>
#include <luabind/exception_handler.hpp>
#include <luabind/shared_ptr_converter.hpp> //< boost::shared_ptr up- and down-casts

#include <halmd/utility/lua_wrapper/any_converter.hpp>
#include <halmd/utility/lua_wrapper/array_converter.hpp>
#include <halmd/utility/lua_wrapper/fixed_vector_converter.hpp>
#include <halmd/utility/lua_wrapper/long_long_converter.hpp>
#include <halmd/utility/lua_wrapper/map_converter.hpp>
#include <halmd/utility/lua_wrapper/optional_converter.hpp>
#include <halmd/utility/lua_wrapper/registry.hpp>
#include <halmd/utility/lua_wrapper/ublas.hpp>
#include <halmd/utility/lua_wrapper/vector_converter.hpp>

/**
 * @namespace halmd::lua_wrapper
 *
 * Lua C++ wrapper extensions
 */

#endif /* ! HALMD_UTILITY_LUA_WRAPPER_LUA_WRAPPER_HPP */
