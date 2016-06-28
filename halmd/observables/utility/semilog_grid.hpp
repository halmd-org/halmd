/*
 * Copyright © 2011-2013  Felix Höfling
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

#ifndef HALMD_OBSERVABLES_UTILITY_SEMILOG_GRID_HPP
#define HALMD_OBSERVABLES_UTILITY_SEMILOG_GRID_HPP

#include <lua.hpp>
#include <vector>

namespace halmd {
namespace observables {
namespace utility {

/**
 * construct a semi-logarithmically spaced grid
 */

class semilog_grid
{
public:
    static void luaopen(lua_State* L);

    // set up grid upon construction, disable decimation by default
    semilog_grid(double start, double stop, unsigned int decimation=0);

    //! returns list of grid points
    std::vector<double> const& value() const
    {
        return grid_;
    }

protected:
    /** list of grid points */
    std::vector<double> grid_;
};

} // namespace observables
} // namespace utility
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_UTILITY_SEMILOG_GRID_HPP */
