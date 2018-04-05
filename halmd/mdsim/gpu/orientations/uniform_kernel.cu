/*
 * Copyright © 2017-2018  Jake Atwell
 * Copyright © 2016       Manuel Dibak
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/gpu/orientations/uniform_kernel.hpp>
#include <halmd/mdsim/positions/lattice_primitive.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/random/gpu/random_number_generator.cuh>
#include <halmd/utility/gpu/thread.cuh>

#include <cmath>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace orientations {
namespace uniform_kernel {

// does rng need to be passed as a reference?
// how is the state stored for other threads?
template <typename float_type, typename rng_type>
__device__ void initialize_orientation(
        fixed_vector<float_type, 2> & u
      , rng_type & rng
      )
{
        //read random number generator state from global device memory
        typename rng_type::state_type state = rng[GTID];

        // pick a random number between 0 and 2π
        float theta =  6.283185 * random::gpu::uniform(rng, state);

        /* select random point on unit circle */
        u[0] = cos(theta);
        u[1] = sin(theta);

        rng[GTID] = state;
}

template <typename float_type, typename rng_type>
__device__ void initialize_orientation(
        fixed_vector<float_type, 3> & u
      , rng_type & rng
      )
{
        //read random number generator state from global device memory
        typename rng_type::state_type state = rng[GTID];

        float theta = random::gpu::uniform(rng, state);
        float phi = random::gpu::uniform(rng, state);
        float pi = 4*atanf(1);

        theta =  acosf(2*theta - 1);
        phi = 2*pi*phi;

        /* select random point on unit sphere */
        u[0] = sin(theta) * cos(phi);
        u[1] = sin(theta) * sin(phi);
        u[2] = cos(theta);



        rng[GTID] = state;
}

template <typename vector_type, typename rng_type>
__global__ void uniform(
    float4* g_u
  , unsigned int npart
  , rng_type rng
)
{
    enum { dimension = vector_type::static_size };
    unsigned int const threads = GTDIM;

    for (unsigned int i = GTID; i < npart; i += threads) {
        // load particle orientation
        vector_type u;
        unsigned int nothing;
#ifdef USE_VERLET_DSFUN
        tie(u, nothing) <<= tie(g_u[i], g_u[i + threads]);
#else
        tie(u, nothing) <<= g_u[i];
#endif

        initialize_orientation(u, rng);

#ifdef USE_VERLET_DSFUN
        tie(g_u[i], g_u[i + threads]) <<= tie(u, nothing);
#else
        g_u[i] <<= tie(u, nothing);
#endif
    }
}

} // namespace uniform_kernel

template <typename rng_type, int dimension>
uniform_wrapper<rng_type, dimension> const uniform_wrapper<rng_type, dimension>::kernel = {
#ifdef USE_VERLET_DSFUN
    uniform_kernel::uniform<fixed_vector<dsfloat, dimension>, rng_type>
#else
    uniform_kernel::uniform<fixed_vector<float, dimension>, rng_type>
#endif
};

//template class lattice_wrapper<close_packed_lattice<fixed_vector<float, 2>, fixed_vector<unsigned int, 2> > >;

template class uniform_wrapper<random::gpu::rand48_rng, 2>;
template class uniform_wrapper<random::gpu::rand48_rng, 3>;
template class uniform_wrapper<random::gpu::mrg32k3a_rng, 2>;
template class uniform_wrapper<random::gpu::mrg32k3a_rng, 3>;
} // namespace mdsim
} // namespace gpu
} // namespace orientations
} // namespace halmd
