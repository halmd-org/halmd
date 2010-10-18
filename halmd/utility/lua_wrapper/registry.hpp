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

#include <boost/function.hpp>
#include <boost/lexical_cast.hpp> //< Lua class names from template parameters
#include <map>
#include <utility>

#include <lua.hpp>

namespace halmd
{
namespace lua_wrapper
{
namespace detail
{

struct registry
{
    typedef boost::function<void (lua_State*)> function_type;
    typedef std::multimap<int, function_type>::const_iterator const_iterator;

    /**
     * Returns a pointer to the singleton registry.
     *
     * Read the following to learn why we use a static local variable:
     *
     * What's the "static initialization order fiasco"?
     * http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12
     */
    static std::multimap<int, function_type>& get()
    {
        static std::multimap<int, function_type> registry; //< ordered map
        return registry;
    }
};

} // namespace detail

/**
   Register Luabind C++ wrapper.

   HALMD stores a singleton registry with functions that wrap C++ classes
   and functions for use with Lua, where every C++ module registers itself
   before program startup. This avoids explicit dependencies between module
   source files and thus eases maintenance.

   register_ receives an integer priority to order registration functions
   according to their integer priority. This is used to register base
   classes before their derived classes, where the priority equals the
   distance of a derived class to its base class.

   For example, to register two classes A and B

   @code
namespace halmd
{

class A {};
class B : A {};

   @endcode

   invoke register_ with priority 0 for class A, and priority 1 for class B

   @code
static void register_class_A(lua_State* L)
{
    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            class_<A, boost::shared_ptr<A> >("A")
        ]
    ];
}

static void register_class_B(lua_State* L)
{
    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            class_<B, boost::shared_ptr<A>, bases<A> >("B")
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0)
    [
        bind(&register_class_A, _1)
    ];
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        bind(&register_class_B, _1)
    ];
}

} // namespace halmd
   @endcode
 */
struct register_
{
    register_(int priority_) : priority_(priority_) {}

    register_ const& operator[](detail::registry::function_type const& function) const
    {
        detail::registry::get().insert(std::make_pair(priority_, function));
        return *this;
    }

private:
    int priority_;
};

/**
 * Register C++ classes and functions with Lua state.
 *
 * This function is called by halmd::script to process the registered Lua
 * registration functions, ordered by the integer priority passed to
 * register_. This binds the C++ classes and functions as Lua userdata.
 *
 * @param L Lua state
 */
inline void open(lua_State* L)
{
    detail::registry::const_iterator it, end = detail::registry::get().end();
    for (it = detail::registry::get().begin(); it != end; ++it) {
        it->second(L);
    }
}

} // namespace lua_wrapper

} // namespace halmd

#endif /* ! HALMD_UTILITY_LUA_WRAPPER_REGISTRY_HPP */
