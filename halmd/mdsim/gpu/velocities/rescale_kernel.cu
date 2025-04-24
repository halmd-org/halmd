/*
 * Copyright © 2025 Giorgia Marcelli
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

#include <halmd/mdsim/gpu/velocities/rescale_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocities {
namespace rescale_kernel {

/**
 * rescale velocities to match target value of total energy, separately for each particle
 */
template <
    typename ptr_type,
    typename vector_type
>
__global__ void rescale(ptr_type g_v, float const* g_en_pot, uint npart, float target_energy)
{
    typedef typename vector_type::value_type float_type;

    for (uint i = GTID; i < npart; i += GTDIM) {
        // read velocity, mass, and potential energy from global memory
        vector_type v;
        float mass;
        tie(v, mass) <<= g_v[i];
        float en_pot = g_en_pot[i];

        // kinetic energy of this particle
        float_type en_kin = mass * inner_prod(v, v) / 2;

        // Compute velocity scaling factor to match target total energy
        float_type scale = sqrtf((target_energy - en_pot) / en_kin);
        v *= scale;   // rescale to match target kinetic energy

        // write back rescaled velocities to global memory
        g_v[i] <<= tie(v, mass);
    }
}

} // namespace rescale_kernel

// Wrapper instantiation
template <int dimension, typename float_type>
rescale_wrapper<dimension, float_type> rescale_wrapper<dimension, float_type>::kernel = {
    rescale_kernel::rescale<ptr_type, fixed_vector<float_type, dimension>>
};

// Explicit instantiations
#ifdef USE_GPU_SINGLE_PRECISION
template class rescale_wrapper<3, float>;
template class rescale_wrapper<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class rescale_wrapper<3, dsfloat>;
template class rescale_wrapper<2, dsfloat>;
#endif

} // namespace velocities
} // namespace gpu
} // namespace mdsim
} // namespace halmd
