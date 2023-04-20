/*
 * Copyright © 2017-2018  Jake Atwell
 * Copyright © 2016       Manuel Dibak
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
 * Copyright © 2023       Jaslo Ziska
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

#include <cmath>

#include <halmd/mdsim/gpu/orientations/uniform_kernel.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/random/gpu/random_number_generator.cuh>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace orientations {
namespace uniform_kernel {

template <typename ptr_type, typename vector_type, typename rng_type>
__global__ void uniform(
    ptr_type g_u
  , unsigned int npart
  , rng_type rng
)
{
    enum { dimension = vector_type::static_size };

    typename rng_type::state_type state = rng[GTID];

    for (unsigned int i = GTID; i < npart; i += GTDIM) {
        // load particle orientation and mass
        vector_type u;
        float mass;
        tie(u, mass) <<= g_u[i];

        if (dimension == 2) {
            // pick a random number between 0 and 2π
            float theta = 2 * M_PI * uniform(rng, state);

            // select random point on unit circle
            u[0] = cos(theta);
            u[1] = sin(theta);
        } else {
            float theta = acosf(2 * uniform(rng, state) - 1);
            float phi = 2 * M_PI * uniform(rng, state);

            u[0] = sin(theta) * cos(phi);
            u[1] = sin(theta) * sin(phi);
            u[2] = cos(theta);
        }

        // store particle orientation and mass
        g_u[i] <<= tie(u, mass);
    }

    rng[GTID] = state;
}

} // namespace uniform_kernel

template <int dimension, typename float_type, typename rng_type>
uniform_wrapper<dimension, float_type, rng_type> const uniform_wrapper<dimension, float_type, rng_type>::kernel = {
    uniform_kernel::uniform<ptr_type, fixed_vector<float_type, dimension>, rng_type>
};

#ifdef USE_GPU_SINGLE_PRECISION
template class uniform_wrapper<3, float, random::gpu::mrg32k3a_rng>;
template class uniform_wrapper<2, float, random::gpu::mrg32k3a_rng>;
template class uniform_wrapper<3, float, random::gpu::rand48_rng>;
template class uniform_wrapper<2, float, random::gpu::rand48_rng>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class uniform_wrapper<3, dsfloat, random::gpu::mrg32k3a_rng, 2>;
template class uniform_wrapper<2, dsfloat, random::gpu::mrg32k3a_rng, 3>;
template class uniform_wrapper<3, dsfloat, random::gpu::rand48_rng>;
template class uniform_wrapper<2, dsfloat, random::gpu::rand48_rng>;
#endif

} // namespace mdsim
} // namespace gpu
} // namespace orientations
} // namespace halmd
