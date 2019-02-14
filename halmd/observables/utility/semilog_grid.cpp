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

#include <functional>

#include <halmd/io/logger.hpp>
#include <halmd/observables/utility/semilog_grid.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace utility {

semilog_grid::semilog_grid(
    double start
  , double stop
  , unsigned int decimation
)
{
    LOG_DEBUG("construct semi-logarithmically spaced grid");
    LOG_DEBUG("start: " << start << ", stop: " << stop << ", decimation: " << decimation);

    // set up semi-linearly spaced grid
    double h = start;
    unsigned int i = 0;
    for (double x = start; x < stop; ) {
        grid_.push_back(x);
        x += h;
        // double grid spacing after every 'decimation' number of points
        if (decimation > 0 && ++i % decimation == 0) {
            h *= 2;
        }
    }
}

void semilog_grid::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("utility")
            [
                class_<semilog_grid, std::shared_ptr<semilog_grid> >("semilog_grid")
                    .def(constructor<
                        double
                      , double
                      , unsigned int
                    >())
                    .property("value", &semilog_grid::value)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_utility_semilog_grid(lua_State* L)
{
    semilog_grid::luaopen(L);
    return 0;
}

} // namespace observables
} // namespace utility
} // namespace halmd
