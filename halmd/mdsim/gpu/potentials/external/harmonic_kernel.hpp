/*
 * Copyright © 2014 Felix Höfling
 * Copyright © 2020 Jaslo Ziska
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_HARMONIC_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_HARMONIC_KERNEL_HPP

#include <halmd/utility/tuple.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {
namespace harmonic_kernel {

/**
 * harmonic external potential
 */
template <int dimension>
class harmonic
{
public:
    typedef fixed_vector<float, dimension> vector_type;

    /**
     * Construct harmonic potential.
     */
    HALMD_GPU_ENABLED harmonic(cudaTextureObject_t t_param) : t_param_(t_param) {}

    /**
     * Fetch parameters from texture cache for this particle species
     */
    HALMD_GPU_ENABLED void fetch_param(unsigned int species);

    /**
     * Compute force and potential for interaction.
     *
     * @param r   particle position reduced to periodic box
     * @returns   tuple of force vector @f$ -\nabla U(\vec r) @f$ and potential @f$ U(\vec r) @f$
     */
    HALMD_GPU_ENABLED tuple<vector_type, float> operator()(vector_type const& r) const;

private:
    /** potential parameters for given particle species */
    vector_type offset_;
    float stiffness_;
    cudaTextureObject_t t_param_;
};

template <int dimension>
HALMD_GPU_ENABLED tuple<typename harmonic<dimension>::vector_type, float>
harmonic<dimension>::operator()(vector_type const& r) const
{
    vector_type dr = r - offset_;
    vector_type force = stiffness_ * dr;
    float en_pot = inner_prod(force, dr) / 2;

    return make_tuple(force, en_pot);
}

} // namespace harmonic_kernel

struct harmonic_wrapper {};

} // namespace external
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_HARMONIC_KERNEL_HPP */
