/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_GPU_DYNAMICS_MEAN_QUARTIC_DISPLACEMENT_HPP
#define HALMD_OBSERVABLES_GPU_DYNAMICS_MEAN_QUARTIC_DISPLACEMENT_HPP

#include <lua.hpp>

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/observables/dynamics/mean_quartic_displacement.hpp>
#include <halmd/observables/gpu/dynamics/tagged_particle.hpp>
#include <halmd/observables/gpu/samples/phase_space.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

/**
 * Mean-quartic displacement
 */
template <int dimension, typename float_type>
class mean_quartic_displacement
{
public:
    typedef gpu::samples::phase_space<dimension, float_type> sample_type;
    typedef double result_type;
    typedef accumulator<result_type> accumulator_type;

    struct defaults
    {
        static unsigned int blocks();
        static unsigned int threads();
    };

    static char const* module_name() { return "mean_quartic_displacement"; }

    static void luaopen(lua_State* L);

    mean_quartic_displacement(
        unsigned int blocks = defaults::blocks()
      , unsigned int threads = defaults::threads()
    );

    /**
     * Compute mean-quartic displacement from two phase space samples
     *
     * @param first  phase space sample at initial time t1
     * @param second phase space sample at later time t2
     * @returns MSD at lag time t2 - t1, averaged over all particles
     */
    accumulator_type compute(sample_type const& first, sample_type const& second);

private:
    typedef observables::dynamics::mean_quartic_displacement<dimension, float> correlate_function_type;

    /** functor for compuation of mean-quartic displacement */
    reduction<tagged_particle<correlate_function_type, dsfloat> > compute_mqd_;
};

} // namespace dynamics
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_DYNAMICS_MEAN_QUARTIC_DISPLACEMENT_HPP */
