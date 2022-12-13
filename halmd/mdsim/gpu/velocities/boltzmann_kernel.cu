/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/algorithm/gpu/reduction.cuh>
#include <halmd/mdsim/gpu/velocities/boltzmann_kernel.hpp>
#include <halmd/random/gpu/normal_distribution.cuh>
#include <halmd/random/gpu/random_number_generator.cuh>
#include <halmd/utility/gpu/thread.cuh>

using namespace halmd::algorithm::gpu;

//
// Maxwell-Boltzmann distribution at accurate temperature
//

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocities {
namespace boltzmann_kernel {

/**
 * generate Maxwell-Boltzmann distributed velocities and reduce velocity
 */
template <
    typename ptr_type
  , typename vector_type
  , typename rng_type
  , typename T
>
__global__ void gaussian(
    ptr_type g_v
  , unsigned int npart
  , unsigned int nplace
  , float temp
  , dsfloat_ptr<T> g_mv
  , dsfloat* g_mv2
  , dsfloat* g_m
  , rng_type rng
)
{
    enum { dimension = vector_type::static_size };
    typedef typename vector_type::value_type float_type;

    fixed_vector<dsfloat, dimension> mv = 0;
    dsfloat mv2 = 0;
    dsfloat m = 0;

    // read random number generator state from global device memory
    typename rng_type::state_type state = rng[GTID];

    // normal distribution parameters
    float const mean = 0.f;
    float const sigma = sqrtf(temp);

    // cache second normal variate for uneven dimensions
    bool cached = false;
    typename vector_type::value_type cache;

    for (uint i = GTID; i < npart; i += GTDIM) {
        vector_type v;
        float mass;
        tie(v, mass) <<= g_v[i];
        for (uint j = 0; j < dimension - 1; j += 2) {
            tie(v[j], v[j + 1]) = normal(rng, state, mean, sigma);
        }
        if (dimension % 2) {
           if ((cached = !cached)) {
               tie(v[dimension - 1], cache) = normal(rng, state, mean, sigma);
           }
           else {
               v[dimension - 1] = cache;
           }
        }
        v /= sqrtf(mass);
        mv += mass * v;
        mv2 += mass * inner_prod(v, v);
        m += mass;
        g_v[i] <<= tie(v, mass);
    }

    // store random number generator state in global device memory
    rng[GTID] = state;

    // compute reduced value for all threads in block
    reduce<ternary_sum_>(mv, mv2, m);

    if (TID < 1) {
        // store block reduced value in global memory
        g_mv[blockIdx.x] = split(mv);
        g_mv2[blockIdx.x] = mv2;
        g_m[blockIdx.x] = m;
    }
}

template <
    typename ptr_type
  , typename vector_type
  , typename T
>
__global__ void shift(ptr_type g_v, uint npart, uint nplace, dsfloat temp, dsfloat_ptr<T> const g_mv, dsfloat const* g_mv2, dsfloat const* g_m, uint size)
{
    enum { dimension = vector_type::static_size };
    typedef typename vector_type::value_type float_type;

    extern __shared__ char __s_array[];
    fixed_vector<dsfloat, dimension>* const s_mv = reinterpret_cast<fixed_vector<dsfloat, dimension>*>(__s_array);
    dsfloat* const s_mv2 = reinterpret_cast<dsfloat*>(&s_mv[size]);
    dsfloat* const s_m = reinterpret_cast<dsfloat*>(&s_mv2[size]);

    fixed_vector<dsfloat, dimension> mv = 0;
    dsfloat mv2 = 0;
    dsfloat m = 0;

    for (uint i = TID; i < size; i += TDIM) {
        s_mv[i] = fixed_vector<dsfloat, dimension>(get<0>(g_mv[i]), get<1>(g_mv[i]));
        s_mv2[i] = g_mv2[i];
        s_m[i] = g_m[i];
    }
    __syncthreads();
    for (uint i = 0; i < size; ++i) {
        mv += s_mv[i];
        mv2 += s_mv2[i];
        m += s_m[i];
    }

    vector_type vcm = vector_type(mv / m);

    for (uint i = GTID; i < npart; i += GTDIM) {
        vector_type v;
        float mass;
        tie(v, mass) <<= g_v[i];
        v -= vcm;
        g_v[i] <<= tie(v, mass);
    }
}

} // namespace boltzmann_kernel

template <int dimension, typename float_type, typename rng_type>
boltzmann_wrapper<dimension, float_type, rng_type> boltzmann_wrapper<dimension, float_type, rng_type>::kernel = {
    boltzmann_kernel::gaussian<ptr_type, fixed_vector<float_type, dimension>, rng_type>
  , boltzmann_kernel::shift<ptr_type, fixed_vector<float_type, dimension> >
};

#ifdef USE_GPU_SINGLE_PRECISION
template class boltzmann_wrapper<3, float, random::gpu::mrg32k3a_rng>;
template class boltzmann_wrapper<2, float, random::gpu::mrg32k3a_rng>;
template class boltzmann_wrapper<3, float, random::gpu::rand48_rng>;
template class boltzmann_wrapper<2, float, random::gpu::rand48_rng>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class boltzmann_wrapper<3, dsfloat, random::gpu::mrg32k3a_rng>;
template class boltzmann_wrapper<2, dsfloat, random::gpu::mrg32k3a_rng>;
template class boltzmann_wrapper<3, dsfloat, random::gpu::rand48_rng>;
template class boltzmann_wrapper<2, dsfloat, random::gpu::rand48_rng>;
#endif

} // namespace mdsim
} // namespace gpu
} // namespace velocities
} // namespace halmd
