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

#ifndef HALMD_OBSERVABLES_BINNED_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_BINNED_PHASE_SPACE_HPP

#include <boost/multi_array.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/utility/signal.hpp>

namespace halmd {
namespace observables {

/**
 * bin phase space sample into spatial cells
 *
 * compute index cell lists from a given phase space sample
 */
template <int dimension>
class binned_phase_space
{
protected:
    typedef signal<void ()> signal_type;

public:
    typedef signal_type::slot_function_type slot_function_type;
    typedef boost::multi_array<double, dimension + 1> position_grid_type;

    static void luaopen(lua_State* L);

    virtual void acquire() = 0;
    virtual connection on_acquire(slot_function_type const& slot) = 0;

    /** returns midpoint positions of cell grid along specified axis */
    virtual std::vector<double> const& position(unsigned int i) const = 0;

    /** returns full grid of cell positions
     *
     *  The last dimension refers to the coordinates of a position vector.
     *
     *  The grid is constructed on function invocation and returned by value.
     */
    virtual position_grid_type position() const;
};

} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_BINNED_PHASE_SPACE_HPP */
