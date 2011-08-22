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

#ifndef HALMD_OBSERVABLES_GPU_DYNAMICS_MEAN_SQUARE_DISPLACEMENT_HPP
#define HALMD_OBSERVABLES_GPU_DYNAMICS_MEAN_SQUARE_DISPLACEMENT_HPP

#include <lua.hpp>

#include <halmd/observables/gpu/dynamics/mean_square_displacement_kernel.hpp>
#include <halmd/observables/gpu/samples/phase_space.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

/**
 * Mean-square displacement
 */
template <int dimension, typename float_type>
class mean_square_displacement
{
public:
    typedef gpu::samples::phase_space<dimension, float_type> sample_type;
    typedef typename sample_type::sample_vector sample_vector_type;
    typedef typename sample_type::vector_type vector_type;
    typedef double result_type;
    typedef accumulator<result_type> accumulator_type;

    struct defaults
    {
        static unsigned int blocks();
        static unsigned int threads();
    };

    static char const* module_name() { return "mean_square_displacement"; }

    static void luaopen(lua_State* L);

    /**
     * @param type particle type for which the computation is done
     */
    mean_square_displacement(
        size_t type
      , unsigned int blocks = defaults::blocks()
      , unsigned int threads = defaults::threads()
    );

    /**
     * Compute mean-square displacement from two phase space samples
     *
     * @param first  phase space sample at initial time t1
     * @param second phase space sample at later time t2
     * @returns MSD at lag time t2 - t1, averaged over all particles of specified type
     */
    accumulator_type compute(sample_type const& first, sample_type const& second);

private:
    typedef accumulator<dsfloat> gpu_accumulator_type;
    typedef typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    typedef cuda::function<void (coalesced_vector_type const*, coalesced_vector_type const*, unsigned int, gpu_accumulator_type*)> compute_function_type;

    static compute_function_type select_compute(unsigned int threads);

    size_t type_;
    unsigned int blocks_;
    unsigned int threads_;
    compute_function_type compute_;
    cuda::vector<gpu_accumulator_type> g_acc_;
    cuda::host::vector<gpu_accumulator_type> h_acc_;
};

} // namespace dynamics
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_DYNAMICS_MEAN_SQUARE_DISPLACEMENT_HPP */
