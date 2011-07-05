/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_HOST_VELOCITY_AUTOCORRELATION_HPP
#define HALMD_OBSERVABLES_HOST_VELOCITY_AUTOCORRELATION_HPP

#include <lua.hpp>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>

namespace halmd {
namespace observables {
namespace host {

/**
 * Velocity autocorrelation.
 */
template <int dimension, typename float_type>
class velocity_autocorrelation
{
public:
    typedef host::samples::phase_space<dimension, float_type> phase_space_type;
    typedef typename phase_space_type::vector_type vector_type;
    typedef typename phase_space_type::sample_vector sample_vector;
    typedef accumulator<float_type> result_type;

    static void luaopen(lua_State* L);

    velocity_autocorrelation() {}
    result_type compute(sample_vector const& first, sample_vector const& second);
};

} // namespace observables
} // namespace host
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_VELOCITY_AUTOCORRELATION_HPP */
