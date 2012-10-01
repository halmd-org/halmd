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

#ifndef HALMD_OBSERVABLES_FIELDS_DENSITY_HPP
#define HALMD_OBSERVABLES_FIELDS_DENSITY_HPP

#include <boost/multi_array.hpp>
#include <lua.hpp>

#include <halmd/utility/signal.hpp>

namespace halmd {
namespace observables {
namespace fields {

/**
 * compute number density field from binned phase space sample
 *
 * The result is multi-dimensional array with the number of particles per bin,
 * normalised by the volume of the binning cell.
 */
template <int dimension>
class density
{
protected:
    typedef signal<void ()> signal_type;

public:
    typedef boost::multi_array<double, dimension> result_type;
    typedef signal_type::slot_function_type slot_function_type;

    static void luaopen(lua_State* L);

    virtual void sample() = 0;
    virtual connection on_sample(slot_function_type const& slot) = 0;

    virtual result_type const& value() const = 0;
};

} // namespace fields
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_FIELDS_DENSITY_HPP */
