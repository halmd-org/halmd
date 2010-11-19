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

#ifndef HALMD_UTILITY_LUA_WRAPPER_ERROR_HPP
#define HALMD_UTILITY_LUA_WRAPPER_ERROR_HPP

#include <luabind/luabind.hpp>

namespace halmd
{
namespace lua_wrapper
{

/**
 * set lua_pcall error handler within current scope
 */
class scoped_pcall_callback
{
public:
    /**
     * set lua_pcall error handler
     *
     * @param pcall_callback lua_pcall error handler
     */
    scoped_pcall_callback(luabind::pcall_callback_fun pcall_callback)
    {
        pcall_callback_ = luabind::get_pcall_callback();

        luabind::set_pcall_callback(pcall_callback);
    }

    /**
     * restore lua_pcall error handler of parent scope
     */
    ~scoped_pcall_callback()
    {
        luabind::set_pcall_callback(pcall_callback_);
    }

private:
    //! lua_pcall error handler of parent scope
    luabind::pcall_callback_fun pcall_callback_;
};

} // namespace lua_wrapper

} // namespace halmd

#endif /* ! HALMD_UTILITY_LUA_WRAPPER_ERROR_HPP */
