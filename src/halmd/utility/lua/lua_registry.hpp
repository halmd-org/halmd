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

#ifndef HALMD_UTILITY_LUA_LUA_REGISTRY_HPP
#define HALMD_UTILITY_LUA_LUA_REGISTRY_HPP

#include <boost/mpl/int.hpp>
#include <boost/parameter.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/type_traits/is_same.hpp>
#include <list>

#include <halmd/utility/lua/lua_include.hpp>

namespace halmd
{

namespace detail
{

namespace tag { struct dimension; } // keyword tag type
namespace tag { struct backend; } // keyword tag type

struct backend_base {};

template <typename T>
struct is_backend
  : boost::is_base_and_derived<backend_base, T> {};

template <typename T>
struct is_dimension
  : boost::is_same<boost::mpl::int_<T::value>, T> {};

} // namespace detail

template <int Dimension>
struct dimension_
  : boost::parameter::template_keyword<
        detail::tag::dimension, boost::mpl::int_<Dimension> > {};

template <typename Backend>
struct backend_
  : boost::parameter::template_keyword<
        detail::tag::backend, Backend> {};

struct gpu_ : detail::backend_base {};
struct host_ : detail::backend_base {};

namespace detail
{

typedef boost::parameter::parameters<
    boost::parameter::optional<
        boost::parameter::deduced<detail::tag::backend>
      , detail::is_backend<boost::mpl::_>
    >
  , boost::parameter::optional<
        detail::tag::dimension
      , detail::is_dimension<boost::mpl::_>
    >
> registry_signature;

struct lua_registry_impl_base
{
    /** wrapper function type */
    typedef void (*value_type)(lua_State*);
    /** registry iterator */
    typedef std::list<value_type>::iterator iterator;
    /** registry const iterator */
    typedef std::list<value_type>::const_iterator const_iterator;
};

template <
    typename Dimension
  , typename Backend
>
struct lua_registry_impl
  : lua_registry_impl_base
{
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

} // namespace detail

/**
 * This class stores a static list of C++ wrapper registration functions.
 *
 * To export C++ classes or functions to Lua, we have to register them
 * using the Luabind C++ wrapper API. The wrappers are implemented as free
 * functions, which are pushed into this registry at program startup using
 * static initialization. The registry is then iterated over to initialize
 * the Lua state.
 */
template <
    typename A0 = boost::parameter::void_
  , typename A1 = boost::parameter::void_
>
class lua_registry
  : public detail::lua_registry_impl_base
{
private:
    // Create ArgumentPack
    typedef typename
        detail::registry_signature::bind<A0, A1>::type args;

    // Extract first logical parameter.
    typedef typename boost::parameter::binding<
      args, detail::tag::dimension, void>::type dimension;

    typedef typename boost::parameter::binding<
      args, detail::tag::backend, void>::type backend;

public:
    /**
     * Returns a pointer to the singleton registry.
     */
    static std::list<value_type>* get()
    {
        return detail::lua_registry_impl<dimension, backend>::get();
    }
};

} // namespace halmd

#endif /* ! HALMD_UTILITY_LUA_LUA_REGISTRY_HPP */
