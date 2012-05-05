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

#include <boost/function.hpp>
#include <luabind/luabind.hpp>
#include <stdexcept>

/**
 * Luabind converter for Boost.Function
 */

namespace halmd { namespace detail {

template <typename T>
class lua_to_cpp_function;

#ifndef HALMD_NO_CXX11

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

#else /* HALMD_NO_CXX11 */

template <typename T0>
class lua_to_cpp_function<T0 ()>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    T0 operator()() const
    {
        try {
            return luabind::call_function<T0>(f_);
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

template <typename T0, typename T1>
class lua_to_cpp_function<T0 (T1)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    T0 operator()(T1 arg1) const
    {
        try {
            return luabind::call_function<T0>(f_, arg1);
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

template <typename T0, typename T1, typename T2>
class lua_to_cpp_function<T0 (T1, T2)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    T0 operator()(T1 arg1, T2 arg2) const
    {
        try {
            return luabind::call_function<T0>(f_, arg1, arg2);
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

template <typename T0, typename T1, typename T2, typename T3>
class lua_to_cpp_function<T0 (T1, T2, T3)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    T0 operator()(T1 arg1, T2 arg2, T3 arg3) const
    {
        try {
            return luabind::call_function<T0>(f_, arg1, arg2, arg3);
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

template <typename T0, typename T1, typename T2, typename T3, typename T4>
class lua_to_cpp_function<T0 (T1, T2, T3, T4)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4) const
    {
        try {
            return luabind::call_function<T0>(f_, arg1, arg2, arg3, arg4);
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

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
class lua_to_cpp_function<T0 (T1, T2, T3, T4, T5)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5) const
    {
        try {
            return luabind::call_function<T0>(f_, arg1, arg2, arg3, arg4, arg5);
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

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
class lua_to_cpp_function<T0 (T1, T2, T3, T4, T5, T6)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6) const
    {
        try {
            return luabind::call_function<T0>(f_, arg1, arg2, arg3, arg4, arg5, arg6);
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

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
class lua_to_cpp_function<T0 (T1, T2, T3, T4, T5, T6, T7)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7) const
    {
        try {
            return luabind::call_function<T0>(f_, arg1, arg2, arg3, arg4, arg5, arg6, arg7);
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

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
class lua_to_cpp_function<T0 (T1, T2, T3, T4, T5, T6, T7, T8)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7, T8 arg8) const
    {
        try {
            return luabind::call_function<T0>(f_, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
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

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
class lua_to_cpp_function<T0 (T1, T2, T3, T4, T5, T6, T7, T8, T9)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7, T8 arg8, T9 arg9) const
    {
        try {
            return luabind::call_function<T0>(f_, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
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

template <>
class lua_to_cpp_function<void ()>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()() const
    {
        try {
            luabind::call_function<void>(f_);
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

template <typename T1>
class lua_to_cpp_function<void (T1)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()(T1 arg1) const
    {
        try {
            luabind::call_function<void>(f_, arg1);
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

template <typename T1, typename T2>
class lua_to_cpp_function<void (T1, T2)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()(T1 arg1, T2 arg2) const
    {
        try {
            luabind::call_function<void>(f_, arg1, arg2);
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

template <typename T1, typename T2, typename T3>
class lua_to_cpp_function<void (T1, T2, T3)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()(T1 arg1, T2 arg2, T3 arg3) const
    {
        try {
            luabind::call_function<void>(f_, arg1, arg2, arg3);
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

template <typename T1, typename T2, typename T3, typename T4>
class lua_to_cpp_function<void (T1, T2, T3, T4)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4) const
    {
        try {
            luabind::call_function<void>(f_, arg1, arg2, arg3, arg4);
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

template <typename T1, typename T2, typename T3, typename T4, typename T5>
class lua_to_cpp_function<void (T1, T2, T3, T4, T5)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5) const
    {
        try {
            luabind::call_function<void>(f_, arg1, arg2, arg3, arg4, arg5);
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

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
class lua_to_cpp_function<void (T1, T2, T3, T4, T5, T6)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6) const
    {
        try {
            luabind::call_function<void>(f_, arg1, arg2, arg3, arg4, arg5, arg6);
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

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
class lua_to_cpp_function<void (T1, T2, T3, T4, T5, T6, T7)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7) const
    {
        try {
            luabind::call_function<void>(f_, arg1, arg2, arg3, arg4, arg5, arg6, arg7);
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

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
class lua_to_cpp_function<void (T1, T2, T3, T4, T5, T6, T7, T8)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7, T8 arg8) const
    {
        try {
            luabind::call_function<void>(f_, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
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

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
class lua_to_cpp_function<void (T1, T2, T3, T4, T5, T6, T7, T8, T9)>
{
public:
    explicit lua_to_cpp_function(luabind::object const& function) : f_(function) {}

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7, T8 arg8, T9 arg9) const
    {
        try {
            luabind::call_function<void>(f_, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
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

#endif /* HALMD_NO_CXX11 */

template <typename Base>
struct lua_function_converter
  : Base
{
public:
    //! compute Lua to C++ conversion score
    template <typename T>
    int match(lua_State* L, luabind::detail::by_value<boost::function<T> > t, int index)
    {
        return lua_isfunction(L, index) ? 0 : Base::match(L, t, index);
    }

    //! compute Lua to C++ conversion score
    template <typename T>
    int match(lua_State* L, luabind::detail::by_const_reference<boost::function<T> > t, int index)
    {
        return lua_isfunction(L, index) ? 0 : Base::match(L, t, index);
    }

    //! convert from Lua to C++
    template <typename T>
    boost::function<T> apply(lua_State* L, luabind::detail::by_value<boost::function<T> > t, int index)
    {
        return lua_isfunction(L, index) ? lua_to_cpp_function<T>(luabind::object(luabind::from_stack(L, index))) : Base::apply(L, t, index);
    }

    //! convert from Lua to C++
    template <typename T>
    boost::function<T> apply(lua_State* L, luabind::detail::by_const_reference<boost::function<T> > t, int index)
    {
        return lua_isfunction(L, index) ? lua_to_cpp_function<T>(luabind::object(luabind::from_stack(L, index))) : Base::apply(L, t, index);
    }

    //! convert from C++ to Lua
    template <typename T>
    void apply(lua_State* L, boost::function<T> const& value)
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
    int match(lua_State* L, luabind::detail::by_value<boost::function<T> > t, int index)
    {
        return Base::match(L, t, index);
    }

    //! compute Lua to C++ conversion score
    template <typename T>
    int match(lua_State* L, luabind::detail::by_const_reference<boost::function<T> > t, int index)
    {
        return Base::match(L, t, index);
    }

    //! convert from Lua to C++
    template <typename T>
    boost::function<T> apply(lua_State* L, luabind::detail::by_value<boost::function<T> > t, int index)
    {
        return Base::apply(L, t, index);
    }

    //! convert from Lua to C++
    template <typename T>
    boost::function<T> apply(lua_State* L, luabind::detail::by_const_reference<boost::function<T> > t, int index)
    {
        return Base::apply(L, t, index);
    }

    template <typename T>
    void apply(lua_State* L, boost::function<T> const& value)
    {
        if (!luabind::detail::class_registry::get_registry(L)->find_class(typeid(boost::function<T>))) {
            luabind::module(L)
            [
                luabind::class_<boost::function<T> >(typeid(boost::function<T>).name())
                    .def("__call", &boost::function<T>::operator())
            ];
        }
        Base::apply(L, value);
    }
};

}} // namespace halmd::detail

namespace luabind {

#ifndef HALMD_NO_CXX11

template <typename R, typename... Args>
struct default_converter<boost::function<R (Args...)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename R, typename... Args>
struct default_converter<boost::function<R (Args...)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename R, typename... Args>
struct default_converter<boost::function<R& (Args...)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename R, typename... Args>
struct default_converter<boost::function<R& (Args...)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

#else /* HALMD_NO_CXX11 */

template <typename T0>
struct default_converter<boost::function<T0 ()> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename T0, typename T1>
struct default_converter<boost::function<T0 (T1)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename T0, typename T1, typename T2>
struct default_converter<boost::function<T0 (T1, T2)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename T0, typename T1, typename T2, typename T3>
struct default_converter<boost::function<T0 (T1, T2, T3)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4>
struct default_converter<boost::function<T0 (T1, T2, T3, T4)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
struct default_converter<boost::function<T0 (T1, T2, T3, T4, T5)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
struct default_converter<boost::function<T0 (T1, T2, T3, T4, T5, T6)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
struct default_converter<boost::function<T0 (T1, T2, T3, T4, T5, T6, T7)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
struct default_converter<boost::function<T0 (T1, T2, T3, T4, T5, T6, T7, T8)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
struct default_converter<boost::function<T0 (T1, T2, T3, T4, T5, T6, T7, T8, T9)> >
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::value_converter> > {};


template <typename T0>
struct default_converter<boost::function<T0 ()> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename T0, typename T1>
struct default_converter<boost::function<T0 (T1)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename T0, typename T1, typename T2>
struct default_converter<boost::function<T0 (T1, T2)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename T0, typename T1, typename T2, typename T3>
struct default_converter<boost::function<T0 (T1, T2, T3)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4>
struct default_converter<boost::function<T0 (T1, T2, T3, T4)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
struct default_converter<boost::function<T0 (T1, T2, T3, T4, T5)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
struct default_converter<boost::function<T0 (T1, T2, T3, T4, T5, T6)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
struct default_converter<boost::function<T0 (T1, T2, T3, T4, T5, T6, T7)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
struct default_converter<boost::function<T0 (T1, T2, T3, T4, T5, T6, T7, T8)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
struct default_converter<boost::function<T0 (T1, T2, T3, T4, T5, T6, T7, T8, T9)> const&>
  : halmd::detail::lua_function_converter<halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> > {};


template <typename T0>
struct default_converter<boost::function<T0& ()> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename T0, typename T1>
struct default_converter<boost::function<T0& (T1)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename T0, typename T1, typename T2>
struct default_converter<boost::function<T0& (T1, T2)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename T0, typename T1, typename T2, typename T3>
struct default_converter<boost::function<T0& (T1, T2, T3)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4>
struct default_converter<boost::function<T0& (T1, T2, T3, T4)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
struct default_converter<boost::function<T0& (T1, T2, T3, T4, T5)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
struct default_converter<boost::function<T0& (T1, T2, T3, T4, T5, T6)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
struct default_converter<boost::function<T0& (T1, T2, T3, T4, T5, T6, T7)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
struct default_converter<boost::function<T0& (T1, T2, T3, T4, T5, T6, T7, T8)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
struct default_converter<boost::function<T0& (T1, T2, T3, T4, T5, T6, T7, T8, T9)> >
  : halmd::detail::cpp_function_converter<luabind::detail::value_converter> {};


template <typename T0>
struct default_converter<boost::function<T0& ()> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

template <typename T0, typename T1>
struct default_converter<boost::function<T0& (T1)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

template <typename T0, typename T1, typename T2>
struct default_converter<boost::function<T0& (T1, T2)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

template <typename T0, typename T1, typename T2, typename T3>
struct default_converter<boost::function<T0& (T1, T2, T3)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4>
struct default_converter<boost::function<T0& (T1, T2, T3, T4)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
struct default_converter<boost::function<T0& (T1, T2, T3, T4, T5)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
struct default_converter<boost::function<T0& (T1, T2, T3, T4, T5, T6)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
struct default_converter<boost::function<T0& (T1, T2, T3, T4, T5, T6, T7)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
struct default_converter<boost::function<T0& (T1, T2, T3, T4, T5, T6, T7, T8)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
struct default_converter<boost::function<T0& (T1, T2, T3, T4, T5, T6, T7, T8, T9)> const&>
  : halmd::detail::cpp_function_converter<luabind::detail::const_ref_converter> {};

#endif /* HALMD_NO_CXX11 */

} // namespace luabind

#endif /* ! HALMD_UTILITY_LUA_FUNCTION_HPP */
