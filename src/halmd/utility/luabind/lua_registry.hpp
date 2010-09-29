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

#ifndef HALMD_UTILITY_LUABIND_LUA_REGISTRY_HPP
#define HALMD_UTILITY_LUABIND_LUA_REGISTRY_HPP

#include <list>

#include <luabind/luabind.hpp>

namespace halmd
{

/**
 * This class stores a static list of C++ wrapper registration functions.
 *
 * To export C++ classes or functions to Lua, we have to register them
 * using the Luabind C++ wrapper API. The wrappers are implemented as free
 * functions, which are pushed into this registry at program startup using
 * static initialization. The registry is then iterated over to initialize
 * the Lua state.
 */
struct lua_registry
{
    /** wrapper function type */
    typedef luabind::scope value_type;
    /** registry iterator */
    typedef std::list<value_type>::iterator iterator;
    /** registry const iterator */
    typedef std::list<value_type>::const_iterator const_iterator;

    /**
     * Returns a pointer to the singleton registry.
     *
     * Read the following to learn why we use a static local variable:
     *
     * What's the "static initialization order fiasco"?
     * http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12
     */
    static std::list<value_type>* get()
    {
        static std::list<value_type> list;
        return &list;
    }
};

} // namespace halmd

#endif /* ! HALMD_UTILITY_LUABIND_LUA_REGISTRY_HPP */
