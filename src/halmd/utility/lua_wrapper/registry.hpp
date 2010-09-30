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

#ifndef HALMD_UTILITY_LUA_WRAPPER_REGISTRY_HPP
#define HALMD_UTILITY_LUA_WRAPPER_REGISTRY_HPP

#include <map>
#include <utility>

#include <luabind/luabind.hpp>

namespace halmd
{
namespace lua_wrapper
{
namespace detail
{

struct registry
{
    typedef std::multimap<int, luabind::scope>::const_iterator const_iterator;

    /**
     * Returns a pointer to the singleton registry.
     *
     * Read the following to learn why we use a static local variable:
     *
     * What's the "static initialization order fiasco"?
     * http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12
     */
    static std::multimap<int, luabind::scope>& get()
    {
        static std::multimap<int, luabind::scope> registry; //< ordered map
        return registry;
    }
};

} // namespace detail

/**
   Register Luabind C++ wrapper.

   HALMD stores a singleton registry with Lua wrappers of C++ classes and
   functions, where every C++ module registers itself before program
   startup. This avoids explicit dependencies between module source
   files and thus eases maintainance.

   register_ receives an integer priority and a luabind::scope, which may
   be one or multiple (concatenated with comma) namespaces, classes or free
   functions. Scopes are ordered in the registry according to their integer
   priority. This is used to register base classes before their derived
   classes. The priority equals the distance of a derived class to its
   base class.

   For example, to register two classes A and B

   @code
namespace halmd
{

class A {};
class B : A {};

   @endcode

   invoke register_ with priority 0 for class A, and priority 1 for class B

   @code
static __attribute__((constructor)) void register_lua()
{
    using namespace luabind;
    lua_wrapper::register_(0)
    [
        namespace_("halmd_wrapper")
        [
            class_<A, shared_ptr<A> >("A")
        ]
    ];
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        namespace_("halmd_wrapper")
        [
            class_<B, shared_ptr<A>, bases<A> >("B")
        ]
    ];
}

} // namespace halmd
   @endcode
 */
struct register_
{
    register_(int priority_) : priority_(priority_) {}

    void operator[](luabind::scope const& scope) const
    {
        detail::registry::get().insert(std::make_pair(priority_, scope));
    }

private:
    int priority_;
};

/**
 * Load registered C++ wrappers into Lua state.
 *
 * This function is called by halmd::script to process all registered
 * Luabind scopes, ordered by the integer priority passed to register_.
 * This binds the C++ classes and functions to Lua objects.
 *
 * @param L Lua state
 */
inline void open(lua_State* L)
{
    detail::registry::const_iterator it, end = detail::registry::get().end();
    for (it = detail::registry::get().begin(); it != end; ++it) {
        luabind::module(L)[it->second];
    }
}

} // namespace lua_wrapper

} // namespace halmd

#endif /* ! HALMD_UTILITY_LUA_WRAPPER_REGISTRY_HPP */
