/*
 * Copyright © 2012  Peter Colberg
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

#ifndef HALMD_UTILITY_LUA_FUNCTION_HPP
#define HALMD_UTILITY_LUA_FUNCTION_HPP

#include <functional>
#include <luabind/luabind.hpp>
#include <stdexcept>

/**
 * Luabind converter for Boost.Function
 */

namespace halmd { namespace detail {

template <typename T>
class lua_to_cpp_function;

template <typename R, typename... Args>
class lua_to_cpp_function<R (Args...)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    R operator()(Args... args) const
    {
        try {
            return luabind::call_function<R>(f_, args...);
        }
        catch (luabind::error const& e) {
            std::string error(lua_tostring(e.state(), -1));
            lua_pop(e.state(), 1);
            throw std::runtime_error(error);
        }
    }

private:
    luabind::object f_;
};

template <typename... Args>
class lua_to_cpp_function<void (Args...)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()(Args... args) const
    {
        try {
            luabind::call_function<void>(f_, args...);
        }
        catch (luabind::error const& e) {
            std::string error(lua_tostring(e.state(), -1));
            lua_pop(e.state(), 1);
            throw std::runtime_error(error);
        }
    }

private:
    luabind::object f_;
};

template <typename Base>
struct lua_function_converter
  : Base
{
public:
    //! compute Lua to C++ conversion score
    template <typename T>
    int match(lua_State* L, luabind::detail::by_value<std::function<T> > t, int index)
    {
        return lua_isfunction(L, index) ? 0 : Base::match(L, t, index);
    }

    //! compute Lua to C++ conversion score
    template <typename T>
    int match(lua_State* L, luabind::detail::by_const_reference<std::function<T> > t, int index)
    {
        return lua_isfunction(L, index) ? 0 : Base::match(L, t, index);
    }

    //! convert from Lua to C++
    template <typename T>
    std::function<T> apply(lua_State* L, luabind::detail::by_value<std::function<T> > t, int index)
    {
        return lua_isfunction(L, index) ? lua_to_cpp_function<T>(luabind::object(luabind::from_stack(L, index))) : Base::apply(L, t, index);
    }

    //! convert from Lua to C++
    template <typename T>
    std::function<T> apply(lua_State* L, luabind::detail::by_const_reference<std::function<T> > t, int index)
    {
        return lua_isfunction(L, index) ? lua_to_cpp_function<T>(luabind::object(luabind::from_stack(L, index))) : Base::apply(L, t, index);
    }

    //! convert from C++ to Lua
    template <typename T>
    void apply(lua_State* L, std::function<T> const& value)
    {
        Base::apply(L, value);
    }
};

template <typename Base>
struct cpp_function_converter
  : Base
{
    //! compute Lua to C++ conversion score
    template <typename T>
    int match(lua_State* L, luabind::detail::by_value<std::function<T> > t, int index)
    {
        return Base::match(L, t, index);
    }

    //! compute Lua to C++ conversion score
    template <typename T>
    int match(lua_State* L, luabind::detail::by_const_reference<std::function<T> > t, int index)
    {
        return Base::match(L, t, index);
    }

    //! convert from Lua to C++
    template <typename T>
    std::function<T> apply(lua_State* L, luabind::detail::by_value<std::function<T> > t, int index)
    {
        return Base::apply(L, t, index);
    }

    //! convert from Lua to C++
    template <typename T>
    std::function<T> apply(lua_State* L, luabind::detail::by_const_reference<std::function<T> > t, int index)
    {
        return Base::apply(L, t, index);
    }

    template <typename T>
    void apply(lua_State* L, std::function<T> const& value)
    {
        if (!luabind::detail::class_registry::get_registry(L)->find_class(typeid(std::function<T>))) {
            luabind::module(L)
            [
                luabind::class_<std::function<T> >()
                    .def("__call", &std::function<T>::operator())
            ];
        }
        Base::apply(L, value);
    }
};

}} // namespace halmd::detail

namespace luabind {

template <typename R, typename... Args>
struct default_converter<std::function<R (Args...)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename R, typename... Args>
struct default_converter<std::function<R (Args...)>&&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename R, typename... Args>
struct default_converter<std::function<R (Args...)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename R, typename... Args>
struct default_converter<std::function<R& (Args...)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename R, typename... Args>
struct default_converter<std::function<R& (Args...)>&&>
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename R, typename... Args>
struct default_converter<std::function<R& (Args...)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

} // namespace luabind

#endif /* ! HALMD_UTILITY_LUA_FUNCTION_HPP */
