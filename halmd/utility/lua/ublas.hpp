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

#ifndef HALMD_UTILITY_LUA_UBLAS_HPP
#define HALMD_UTILITY_LUA_UBLAS_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp> // row()
#include <boost/numeric/ublas/vector.hpp>

#include <luabind/luabind.hpp>

#if LUA_VERSION_NUM < 502
# define luaL_len lua_objlen
#endif

namespace luabind {

/**
 * Luabind converter for Boost uBLAS vector
 */
template <typename T, typename A>
struct default_converter<boost::numeric::ublas::vector<T, A> >
  : native_converter_base<boost::numeric::ublas::vector<T, A> >
{
    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        return lua_type(L, index) == LUA_TTABLE ? 0 : -1;
    }

    //! convert from Lua to C++
    boost::numeric::ublas::vector<T, A> from(lua_State* L, int index)
    {
        std::size_t size = luaL_len(L, index);
        boost::numeric::ublas::vector<T, A> v(size);
        object table(from_stack(L, index));
        for (std::size_t i = 0; i < v.size(); ++i) {
            v[i] = object_cast<T>(table[i + 1]);
        }
        return v;
    }

    //! convert from C++ to Lua
    void to(lua_State* L, boost::numeric::ublas::vector<T, A> const& v)
    {
        object table = newtable(L);
        for (std::size_t i = 0; i < v.size(); ++i) {
            // default_converter<T> only invoked with reference wrapper
            table[i + 1] = boost::cref(v[i]);
        }
        table.push(L);
    }
};

template <typename T, typename A>
struct default_converter<boost::numeric::ublas::vector<T, A> const&>
  : default_converter<boost::numeric::ublas::vector<T, A> > {};

/**
 * Luabind converter for Boost uBLAS matrix
 */
template <typename T, typename F, typename A>
struct default_converter<boost::numeric::ublas::matrix<T, F, A> >
  : native_converter_base<boost::numeric::ublas::matrix<T, F, A> >
{
    typedef boost::numeric::ublas::vector<T, A> vector_type;
    typedef boost::numeric::ublas::matrix<T, F, A> matrix_type;

    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        return lua_type(L, index) == LUA_TTABLE ? 0 : -1;
    }

    //! convert from Lua to C++
    matrix_type from(lua_State* L, int index)
    {
        object table(from_stack(L, index));
        std::size_t rows = luaL_len(L, index);
        std::size_t cols = (rows > 0) ? object_cast<vector_type>(table[1]).size() : 0;
        matrix_type m(rows, cols);
        for (std::size_t i = 0; i < m.size1(); ++i) {
            row(m, i) = object_cast<vector_type>(table[i + 1]);
        }
        return m;
    }

    //! convert from C++ to Lua
    void to(lua_State* L, matrix_type const& m)
    {
        object table = newtable(L);
        for (std::size_t i = 0; i < m.size1(); ++i) {
            vector_type r = row(m, i);
            // default_converter<T> only invoked with reference wrapper
            table[i + 1] = boost::cref(r);
        }
        table.push(L);
    }
};

template <typename T, typename F, typename A>
struct default_converter<boost::numeric::ublas::matrix<T, F, A> const&>
  : default_converter<boost::numeric::ublas::matrix<T, F, A> > {};

/**
 * Luabind converter for Boost uBLAS unbounded storage array
 */
template <typename T, typename A>
struct default_converter<boost::numeric::ublas::unbounded_array<T, A> >
  : native_converter_base<boost::numeric::ublas::unbounded_array<T, A> >
{
    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        return lua_type(L, index) == LUA_TTABLE ? 0 : -1;
    }

    //! convert from Lua to C++
    boost::numeric::ublas::unbounded_array<T, A> from(lua_State* L, int index)
    {
        std::size_t size = luaL_len(L, index);
        boost::numeric::ublas::unbounded_array<T, A> v(size);
        object table(from_stack(L, index));
        for (std::size_t i = 0; i < v.size(); ++i) {
            v[i] = object_cast<T>(table[i + 1]);
        }
        return v;
    }

    //! convert from C++ to Lua
    void to(lua_State* L, boost::numeric::ublas::unbounded_array<T, A> const& v)
    {
        object table = newtable(L);
        for (std::size_t i = 0; i < v.size(); ++i) {
            // default_converter<T> only invoked with reference wrapper
            table[i + 1] = boost::cref(v[i]);
        }
        table.push(L);
    }
};

template <typename T, typename A>
struct default_converter<boost::numeric::ublas::unbounded_array<T, A> const&>
  : default_converter<boost::numeric::ublas::unbounded_array<T, A> > {};

} // namespace luabind

#if LUA_VERSION_NUM < 502
# undef luaL_len
#endif

#endif /* ! HALMD_UTILITY_LUA_UBLAS_HPP */
