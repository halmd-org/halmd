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

#ifndef HALMD_OBSERVABLES_HOST_DYNAMICS_MEAN_SQUARE_DISPLACEMENT_HPP
#define HALMD_OBSERVABLES_HOST_DYNAMICS_MEAN_SQUARE_DISPLACEMENT_HPP

#include <lua.hpp>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/dynamics/mean_square_displacement.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace dynamics {

/**
 * Mean-square displacement
 */
template <int dimension, typename float_type>
class mean_square_displacement
{
public:
    typedef host::samples::phase_space<dimension, float_type> sample_type;
    typedef typename sample_type::vector_type vector_type;
    typedef double result_type;
    typedef accumulator<result_type> accumulator_type;

    static void luaopen(lua_State* L);

    /**
     * Compute mean-square displacement from two phase space samples
     *
     * @param first  phase space sample at initial time t1
     * @param second phase space sample at later time t2
     * @returns MSD at lag time t2 - t1, averaged over all particles of specified type
     */
    accumulator_type compute(sample_type const& first, sample_type const& second);

private:
    typedef observables::dynamics::mean_square_displacement<dimension, float_type> correlate_function_type;
};

} // namespace dynamics
} // namespace host
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_DYNAMICS_MEAN_SQUARE_DISPLACEMENT_HPP */
