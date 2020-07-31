/*
 * Copyright © 2008-2014 Felix Höfling
 * Copyright © 2014      Nicolas Höft
 * Copyright © 2008-2011 Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_FORCES_EXTERNAL_KERNEL_CUH
#define HALMD_MDSIM_GPU_FORCES_EXTERNAL_KERNEL_CUH

#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/forces/external_kernel.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {
namespace external_kernel {

/**
 * Compute pair forces, potential energy, and stress tensor for all particles
 */
template <
    bool do_aux               //< compute auxiliary variables in addition to force
  , typename vector_type
  , typename potential_type
  , typename gpu_vector_type
>
__global__ void compute(
    potential_type potential
  , float4 const* g_r
  , gpu_vector_type* g_f
  , float* g_en_pot
  , float* g_stress_pot
  , vector_type box_length
  , bool force_zero
)
{
    enum { dimension = vector_type::static_size };
    typedef typename vector_type::value_type value_type;
    typedef typename type_traits<dimension, float>::stress_tensor_type stress_tensor_type;
    unsigned int const i = GTID;

    // load particle associated with this thread
    unsigned int species;
    vector_type r;
    tie(r, species) <<= g_r[i];

    // enforce periodic boundary conditions
    box_kernel::reduce_periodic(r, box_length);

    // fetch potential
    potential.fetch(species);

    // evaluate force acting on this particle
    // and contribution to potential energy
    vector_type f;
    value_type en_pot;
    tie(f, en_pot) = potential(r);

    // add previous force if not set to zero (i.e., this is not the first contribution)
    if (!force_zero) {
        f += static_cast<vector_type>(g_f[i]);
    }
    // write results to global memory
    g_f[i] = static_cast<vector_type>(f);

    // process auxiliary variables if requested
    if (do_aux) {
        // contribution to stress tensor
        stress_tensor_type stress_pot = 0; // = make_stress_tensor(r, f); FIXME

        // add previous results for auxiliary variables if force is not set to zero
        if (!force_zero) {
            en_pot += g_en_pot[i];
            stress_pot += read_stress_tensor<stress_tensor_type>(g_stress_pot + i, GTDIM);
        }
        // write results to global memory
        g_en_pot[i] = en_pot;
        write_stress_tensor(g_stress_pot + i, stress_pot, GTDIM);
    }

}

} // namespace external_kernel

template <int dimension, typename potential_type>
external_wrapper<dimension, potential_type>
external_wrapper<dimension, potential_type>::kernel = {
    external_kernel::compute<false, fixed_vector<float, dimension>, potential_type>
  , external_kernel::compute<true, fixed_vector<float, dimension>, potential_type>
};

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_EXTERNAL_KERNEL_CUH */
