/*
 * Copyright © 2008-2011  Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_DYNAMICS_VELOCITY_AUTOCORRELATION_KERNEL_HPP
#define HALMD_OBSERVABLES_GPU_DYNAMICS_VELOCITY_AUTOCORRELATION_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

template <int dimension, unsigned int threads>
struct velocity_autocorrelation_wrapper
{
    typedef typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    typedef typename mdsim::type_traits<dimension, float>::vector_type vector_type;
    typedef accumulator<dsfloat> accumulator_type;

    cuda::function<void (coalesced_vector_type const*, coalesced_vector_type const*, unsigned int, accumulator_type*)> compute;
    static velocity_autocorrelation_wrapper const wrapper;
};

} // namespace dynamics
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_DYNAMICS_VELOCITY_AUTOCORRELATION_KERNEL_HPP */
