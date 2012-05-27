/*
 * Copyright © 2012 Peter Colberg
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
 * Luabind converter for std::function
 */

namespace halmd {
namespace detail {

template <typename F, typename R, typename... Args>
class cpp_function_converter
  : public luabind::detail::default_converter_generator<F>::type
{
private:
    typedef typename luabind::detail::default_converter_generator<F>::type _Base;

public:
    /**
     * Convert from Lua to C++.
     */
    template <typename T>
    std::function<R (Args...)> apply(lua_State* L, T t, int index)
    {
        return _Base::apply(L, t, index);
    }

    /**
     * Convert from C++ to Lua.
     */
    void apply(lua_State* L, std::function<R (Args...)> const& value)
    {
        if (!luabind::detail::class_registry::get_registry(L)->find_class(typeid(std::function<R (Args...)>))) {
            luabind::module(L)
            [
                luabind::class_<std::function<R (Args...)> >()
                    .def("__call", &std::function<R (Args...)>::operator())
            ];
        }
        _Base::apply(L, value);
    }
};

template <typename F, typename R, typename... Args>
class lua_function_converter
  : public cpp_function_converter<F, R, Args...>
{
private:
    typedef cpp_function_converter<F, R, Args...> _Base;

public:
    /**
     * Compute Lua to C++ conversion score.
     */
    template <typename T>
    int match(lua_State* L, T t, int index)
    {
        return lua_isfunction(L, index) ? 0 : _Base::match(L, t, index);
    }

    /**
     * Convert from Lua to C++.
     */
    template <typename T>
    std::function<R (Args...)> apply(lua_State* L, T t, int index)
    {
        if (lua_isfunction(L, index)) {
            luabind::object function(luabind::from_stack(L, index));
            return [=](Args... args) -> R {
                try {
                    return luabind::call_function<R>(function, args...);
                }
                catch (luabind::error const& e) {
                    std::string error(lua_tostring(e.state(), -1));
                    lua_pop(e.state(), 1);
                    throw std::runtime_error(error);
                }
            };
        }
        return _Base::apply(L, t, index);
    }

    /**
     * Convert from C++ to Lua.
     */
    template <typename T>
    void apply(lua_State* L, std::function<T> const& value)
    {
        _Base::apply(L, value);
    }
};

template <typename F, typename... Args>
class lua_function_converter<F, void, Args...>
  : public cpp_function_converter<F, void, Args...>
{
private:
    typedef cpp_function_converter<F, void, Args...> _Base;

public:
    /**
     * Compute Lua to C++ conversion score.
     */
    template <typename T>
    int match(lua_State* L, T t, int index)
    {
        return lua_isfunction(L, index) ? 0 : _Base::match(L, t, index);
    }

    /**
     * Convert from Lua to C++.
     */
    template <typename T>
    std::function<void (Args...)> apply(lua_State* L, T t, int index)
    {
        if (lua_isfunction(L, index)) {
            luabind::object function(luabind::from_stack(L, index));
            return [=](Args... args) {
                try {
                    luabind::call_function<void>(function, args...);
                }
                catch (luabind::error const& e) {
                    std::string error(lua_tostring(e.state(), -1));
                    lua_pop(e.state(), 1);
                    throw std::runtime_error(error);
                }
            };
        }
        return _Base::apply(L, t, index);
    }

    /**
     * Convert from C++ to Lua.
     */
    template <typename T>
    void apply(lua_State* L, std::function<T> const& value)
    {
        _Base::apply(L, value);
    }
};

} // namespace detail
} // namespace halmd

namespace luabind {

template <typename R, typename... Args>
struct default_converter<std::function<R (Args...)> >
  : halmd::detail::lua_function_converter<std::function<R (Args...)>, R, Args...> {};

template <typename R, typename... Args>
struct default_converter<std::function<R (Args...)>&&>
  : halmd::detail::lua_function_converter<std::function<R (Args...)>&&, R, Args...> {};

template <typename R, typename... Args>
struct default_converter<std::function<R (Args...)> const&>
  : halmd::detail::lua_function_converter<std::function<R (Args...)> const&, R, Args...> {};

template <typename R, typename... Args>
struct default_converter<std::function<R& (Args...)> >
  : halmd::detail::cpp_function_converter<std::function<R& (Args...)>, R&, Args...> {};

template <typename R, typename... Args>
struct default_converter<std::function<R& (Args...)>&&>
  : halmd::detail::cpp_function_converter<std::function<R& (Args...)>&&, R&, Args...> {};

template <typename R, typename... Args>
struct default_converter<std::function<R& (Args...)> const&>
  : halmd::detail::cpp_function_converter<std::function<R& (Args...)> const&, R&, Args...> {};

} // namespace luabind

#endif /* ! HALMD_UTILITY_LUA_FUNCTION_HPP */
