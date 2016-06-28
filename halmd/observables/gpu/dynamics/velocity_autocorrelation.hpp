/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_GPU_DYNAMICS_VELOCITY_AUTOCORRELATION_HPP
#define HALMD_OBSERVABLES_GPU_DYNAMICS_VELOCITY_AUTOCORRELATION_HPP

#include <lua.hpp>

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/observables/dynamics/velocity_autocorrelation.hpp>
#include <halmd/observables/gpu/dynamics/tagged_particle.hpp>
#include <halmd/observables/gpu/samples/phase_space.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

/**
 * Mean-square displacement
 */
template <int dimension, typename float_type>
class velocity_autocorrelation
{
public:
    typedef gpu::samples::phase_space<dimension, float_type> sample_type;
    typedef double result_type;

    struct defaults
    {
        static unsigned int blocks();
        static unsigned int threads();
    };

    static char const* module_name() { return "velocity_autocorrelation"; }

    static void luaopen(lua_State* L);

    velocity_autocorrelation(
        unsigned int blocks = defaults::blocks()
      , unsigned int threads = defaults::threads()
    );

    /**
     * Compute velocity autocorrelation from two phase space samples
     *
     * @param first  phase space sample at initial time t1
     * @param second phase space sample at later time t2
     * @param result returns MSD at lag time t2 - t1, averaged over all particles
     */
    void operator() (sample_type const& first, sample_type const& second, accumulator<result_type>& result);

private:
    typedef observables::dynamics::velocity_autocorrelation<dimension, float> correlate_function_type;

    /** functor for compuation of mean-quartic displacement */
    reduction<tagged_particle<correlate_function_type, dsfloat> > compute_vacf_;
};

} // namespace dynamics
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_DYNAMICS_VELOCITY_AUTOCORRELATION_HPP */
