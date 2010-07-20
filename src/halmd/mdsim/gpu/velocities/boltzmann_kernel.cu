/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/algorithm/gpu/reduction.cuh>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/mdsim/gpu/velocities/boltzmann_kernel.hpp>
#include <halmd/random/gpu/normal_distribution.cuh>
#include <halmd/random/gpu/random_number_generator.cuh>
#include <halmd/utility/gpu/thread.cuh>

using namespace halmd::algorithm::gpu;
using namespace halmd::numeric::gpu::blas;
using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::random::gpu;

//
// Maxwell-Boltzmann distribution at accurate temperature
//

namespace halmd
{
namespace mdsim { namespace gpu { namespace velocities
{
namespace boltzmann_kernel
{

// random number generator parameters
static __constant__ random_number_generator rng;

/**
 * generate Maxwell-Boltzmann distributed velocities and reduce velocity
 */
template <
    typename vector_type
  , typename rng_type
  , int threads
  , typename T
>
__global__ void gaussian(float4* g_v, uint npart, uint nplace, float temp, T* g_vcm, dsfloat* g_vv)
{
    enum { dimension = vector_type::static_size };

    extern __shared__ char __s_array[]; // CUDA 3.0/3.1 breaks template __shared__ type
    vector<dsfloat, dimension>* const s_vcm = reinterpret_cast<vector_type*>(__s_array);
    dsfloat* const s_vv = reinterpret_cast<dsfloat*>(&s_vcm[TDIM]);

    vector<dsfloat, dimension> vcm = 0;
    dsfloat vv = 0;

    // read random number generator state from global device memory
    typename rng_type::state_type state = get<rng_type>(rng)[GTID];

    // normal distribution parameters
    float const mean = 0.f;
    float const sigma = sqrt(temp);

    // cache second normal variate for uneven dimensions
    bool cached = false;
    typename vector_type::value_type cache;

    for (uint i = GTID; i < npart; i += GTDIM) {
        vector_type v;
        unsigned int tag;
#ifdef USE_VERLET_DSFUN
        tie(v, tag) = untagged<vector_type>(g_v[i], g_v[i + nplace]);
#else
        tie(v, tag) = untagged<vector_type>(g_v[i]);
#endif
        for (uint j = 0; j < dimension - 1; j += 2) {
            tie(v[j], v[j + 1]) = normal(get<rng_type>(rng), state, mean, sigma);
        }
        if (dimension % 2) {
           if ((cached = !cached)) {
               tie(v[dimension - 1], cache) = normal(get<rng_type>(rng), state, mean, sigma);
           }
           else {
               v[dimension - 1] = cache;
           }
        }
        vcm += v;
        vv += inner_prod(v, v);
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + nplace]) = tagged(v, tag);
#else
        g_v[i] = tagged(v, tag);
#endif
    }

    // store random number generator state in global device memory
    get<rng_type>(rng)[GTID] = state;

    // reduced values for this thread
    s_vcm[TID] = vcm;
    s_vv[TID] = vv;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<threads / 2, complex_sum_>(vcm, vv, s_vcm, s_vv);

    if (TID < 1) {
        // store block reduced value in global memory
        tie(g_vcm[blockIdx.x], g_vcm[blockIdx.x + BDIM]) = split(vcm);
        g_vv[blockIdx.x] = vv;
    }
}

template <
    typename vector_type
  , typename T
>
__global__ void shift_rescale(float4* g_v, uint npart, uint nplace, dsfloat temp, T const* g_vcm, dsfloat const* g_vv, uint size)
{
    enum { dimension = vector_type::static_size };
    typedef typename vector_type::value_type float_type;

    extern __shared__ char __s_array[]; // CUDA 3.0/3.1 breaks template __shared__ type
    vector<dsfloat, dimension>* const s_vcm = reinterpret_cast<vector_type*>(__s_array);
    dsfloat* const s_vv = reinterpret_cast<dsfloat*>(&s_vcm[size]);

    vector<dsfloat, dimension> vcm = 0;
    dsfloat vv = 0;

    // compute mean center of mass velocity from block reduced values
    for (uint i = TID; i < size; i += TDIM) {
        s_vcm[i] = vector_type(g_vcm[i], g_vcm[i + size]);
        s_vv[i] = g_vv[i];
    }
    __syncthreads();
    for (uint i = 0; i < size; ++i) {
        vcm += s_vcm[i];
        vv += s_vv[i];
    }
    vcm /= npart;
    vv /= npart;

    // center velocities around origin, then rescale to exactly
    // match the desired temperature;
    // temp = vv / dimension
    // vv changes to vv - v_cm^2 after shifting

    vv -= inner_prod(vcm, vcm);
    float_type coeff = sqrt(temp * static_cast<int>(dimension) / vv);

    for (uint i = GTID; i < npart; i += GTDIM) {
        vector_type v;
        unsigned int tag;
#ifdef USE_VERLET_DSFUN
        tie(v, tag) = untagged<vector_type>(g_v[i], g_v[i + nplace]);
#else
        tie(v, tag) = untagged<vector_type>(g_v[i]);
#endif
        v -= vcm;
        v *= coeff;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + nplace]) = tagged(v, tag);
#else
        g_v[i] = tagged(v, tag);
#endif
    }
}

} // namespace boltzmann_kernel

template <int dimension, typename float_type, typename rng_type>
boltzmann_wrapper<dimension, float_type, rng_type> const boltzmann_wrapper<dimension, float_type, rng_type>::kernel = {
    boltzmann_kernel::gaussian<vector<float_type, dimension>, rng_type, 32>
  , boltzmann_kernel::gaussian<vector<float_type, dimension>, rng_type, 64>
  , boltzmann_kernel::gaussian<vector<float_type, dimension>, rng_type, 128>
  , boltzmann_kernel::gaussian<vector<float_type, dimension>, rng_type, 256>
  , boltzmann_kernel::gaussian<vector<float_type, dimension>, rng_type, 512>
  , boltzmann_kernel::shift_rescale<vector<float_type, dimension> >
  , get<rng_type>(boltzmann_kernel::rng)
};

#ifdef USE_VERLET_DSFUN
template class boltzmann_wrapper<3, dsfloat, random::gpu::rand48_rng>;
template class boltzmann_wrapper<2, dsfloat, random::gpu::rand48_rng>;
#else
template class boltzmann_wrapper<3, float, random::gpu::rand48_rng>;
template class boltzmann_wrapper<2, float, random::gpu::rand48_rng>;
#endif

}}} // namespace mdsim::gpu::velocities

} // namespace halmd
