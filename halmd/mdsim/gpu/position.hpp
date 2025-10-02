/*
 * Copyright © 2016 Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_MDSIM_GPU_POSITION_HPP
#define HALMD_MDSIM_GPU_POSITION_HPP

#include <halmd/utility/cache.hpp>

#include <lua.hpp>

#include <vector>

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * abstract class that defines the interface of position modules
 */
class position
{
public:
    virtual ~position() {}

    /** position lists */
    virtual void set() = 0;

    /** Lua bindings */
    static void luaopen(lua_State* L);
};

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POSITION_HPP */
