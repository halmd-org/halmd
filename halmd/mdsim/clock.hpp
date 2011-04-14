/*
 * Copyright © 2011  Felix Höfling
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

#ifndef HALMD_MDSIM_CLOCK_HPP
#define HALMD_MDSIM_CLOCK_HPP

#include <lua.hpp>
#include <stdint.h>

namespace halmd
{
namespace mdsim
{

/**
 * The module contains the current simulation step and time.
 *
 * The main purpose of this class is to separate it from the
 * dimension-dependent mdsim::core module.
 **/
class clock
{
public:
    static void luaopen(lua_State* L);

    clock();

    //! advance clock by one step and given time increment
    void advance(double timestep)
    {
        ++step_;
        time_ += timestep;
    }

    //! returns MD step counter
    uint64_t step() const
    {
        return step_;
    }

    //! returns MD time
    double time() const
    {
        return time_;
    }

private:
    /** step counter */
    uint64_t step_;
    /** simulation time */
    double time_;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_CLOCK_HPP */
