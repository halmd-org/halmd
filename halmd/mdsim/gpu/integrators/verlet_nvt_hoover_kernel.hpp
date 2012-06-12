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

#ifndef HALMD_MDSIM_GPU_INTEGRATOR_VERLET_NVT_HOOVER_KERNEL_HPP
#define HALMD_MDSIM_GPU_INTEGRATOR_VERLET_NVT_HOOVER_KERNEL_HPP

#include <halmd/config.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type>
struct verlet_nvt_hoover_wrapper
{
    typedef fixed_vector<float, dimension> vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;

    cuda::function <void (float4*, coalesced_vector_type*, float4*, coalesced_vector_type const*, float, float_type, vector_type)> integrate;
    cuda::function <void (float4*, coalesced_vector_type const*, float)> finalize;
    cuda::function <void (float4*, float_type)> rescale;

    static verlet_nvt_hoover_wrapper const kernel;
};

/**
 * Compute total kinetic energy.
 */
template <int dimension, typename float_type>
class kinetic_energy
{
public:
    /** element pointer type of input array */
    typedef float4 const* iterator;

    /**
     * Initialise kinetic energy to zero.
     */
    kinetic_energy() : mv2_(0) {}

    /**
     * Accumulate kinetic energy of a particle.
     */
    HALMD_GPU_ENABLED void operator()(float4 const& velocity)
    {
        fixed_vector<float, dimension> v;
        float mass;
        tie(v, mass) <<= velocity;
        mv2_ += mass * inner_prod(v, v);
    }

    /**
     * Accumulate kinetic energy of another accumulator.
     */
    HALMD_GPU_ENABLED void operator()(kinetic_energy const& acc)
    {
        mv2_ += acc.mv2_;
    }

    /**
     * Returns total kinetic energy.
     */
    float_type operator()() const
    {
        return 0.5 * mv2_;
    }

private:
    /** sum over mass × square of velocity vector */
    float_type mv2_;
};

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATOR_VERLET_NVT_HOOVER_KERNEL_HPP */
