/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_MDSIM_NEIGHBOUR_HPP
#define HALMD_MDSIM_NEIGHBOUR_HPP

#include <lua.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace mdsim {

class neighbour
{
public:
    typedef halmd::signal<void ()> signal_type;
    typedef signal_type::slot_function_type slot_function_type;

    static void luaopen(lua_State* L);

    neighbour() {}
    virtual ~neighbour() {}
    virtual void update() = 0;
    /** add slot to be run before neighbour list update */
    virtual connection on_update(slot_function_type const& slot) = 0;
};

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_NEIGHBOUR_HPP */
