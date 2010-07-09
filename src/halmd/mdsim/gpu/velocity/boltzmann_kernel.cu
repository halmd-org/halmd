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
#include <halmd/mdsim/gpu/velocity/boltzmann_kernel.cuh>
#include <halmd/numeric/gpu/blas/dsfloat.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>
#include <halmd/utility/gpu/thread.cuh>

using namespace halmd::numeric::gpu::blas;
using namespace halmd::mdsim::gpu::particle_kernel;

//
// Maxwell-Boltzmann distribution at accurate temperature
//

namespace halmd
{
namespace mdsim { namespace gpu { namespace velocity
{
namespace boltzmann_kernel
{

enum { BLOCKS = boltzmann_wrapper<>::BLOCKS };
enum { THREADS = boltzmann_wrapper<>::THREADS };

/**
 * generate Maxwell-Boltzmann distributed velocities and reduce velocity
 */
template <typename vector_type, typename T>
__global__ void gaussian(float4* g_v, uint npart, uint nplace, float temp, T* g_vcm)
{
    enum { dimension = vector_type::static_size };
    __shared__ vector_type s_vcm[THREADS];
    vector_type vcm = 0;

//     // read random number generator state from global device memory
//     rand48::state_type state = rand48::g_state[GTID];
//
//     for (uint i = GTID; i < npart; i += GTDIM) {
//         T v;
//         rand48::gaussian(v, temp, state);
//         g_v[i] = v;
// #ifdef USE_VERLET_DSFUN
//         g_v[i + nplace] = vector<float, dimension>(0);
// #endif
//         vcm += vector<float, dimension>(v);
//     }
//     // store random number generator state in global device memory
//     rand48::g_state[GTID] = state;

    // reduced value for this thread
    s_vcm[TID] = vcm;
    __syncthreads();

    // compute reduced value for all threads in block
//     reduce<THREADS / 2, sum_>(vcm, s_vcm);

    if (TID < 1) {
        // store block reduced value in global memory
        g_vcm[blockIdx.x] = static_cast<vector<float, dimension> >(vcm);
#ifdef USE_VERLET_DSFUN
        g_vcm[blockIdx.x + BDIM] = dsfloat_lo(vcm);
#endif
    }
}

/**
 * set center of mass velocity to zero and reduce squared velocity
 */
template <typename vector_type, typename T>
__global__ void shift_velocity(float4* g_v, uint npart, uint nplace, T const* g_vcm, dsfloat* g_vv)
{
    enum { dimension = vector_type::static_size };
    __shared__ vector_type s_vcm[BLOCKS];
    __shared__ dsfloat s_vv[THREADS];
    vector_type vcm = 0;
    dsfloat vv = 0;

    // compute mean center of mass velocity from block reduced values
    for (uint i = TID; i < BLOCKS; i += TDIM) {
#ifdef USE_VERLET_DSFUN
        s_vcm[i] = vector_type(g_vcm[i], g_vcm[i + BDIM]);
#else
        s_vcm[i] = g_vcm[i];
#endif
    }
    __syncthreads();
    for (uint i = 0; i < BLOCKS; ++i) {
        vcm += s_vcm[i];
    }
    vcm /= npart;

    for (uint i = GTID; i < npart; i += GTDIM) {
#ifdef USE_VERLET_DSFUN
        vector_type v; // FIXME (g_v[i], g_v[i + nplace]);
#else
        vector_type v = g_v[i];
#endif
        v -= vcm;
        g_v[i] = tagged(static_cast<vector<float, dimension> >(v), /* FIXME */0);
#ifdef USE_VERLET_DSFUN
        g_v[i + nplace] = tagged(dsfloat_lo(v), /* FIXME */0);
#endif
        vv += inner_prod(v, v);
    }
    // reduced value for this thread
    s_vv[TID] = vv;
    __syncthreads();

    // compute reduced value for all threads in block
//     reduce<THREADS / 2, sum_>(vv, s_vv);

    if (TID < 1) {
        // store block reduced value in global memory
        g_vv[blockIdx.x] = vv;
    }
}

/**
 * rescale velocities to accurate temperature
 */
template <typename vector_type>
__global__ void scale_velocity(float4* g_v, uint npart, uint nplace, dsfloat const* g_vv, dsfloat temp)
{
    enum { dimension = vector_type::static_size };
    typedef typename vector_type::value_type float_type;
    __shared__ dsfloat s_vv[THREADS];
    dsfloat vv = 0;

    // compute squared velocity sum from block reduced values
    for (uint i = TID; i < BLOCKS; i += TDIM) {
        s_vv[i] = g_vv[i];
    }
    __syncthreads();
    for (uint i = 0; i < BLOCKS; ++i) {
        vv += s_vv[i];
    }

    float_type coeff = sqrt(temp * static_cast<int>(dimension) * (static_cast<float_type>(npart) / vv));

    for (uint i = GTID; i < npart; i += GTDIM) {
        vector_type v;
        unsigned int tag;
#ifdef USE_VERLET_DSFUN
        tie(v, tag) = untagged<vector_type>(g_v[i], g_v[i + nplace]);
#else
        tie(v, tag) = untagged<vector_type>(g_v[i]);
#endif
        v *= coeff;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + nplace]) = tagged(v, tag);
#else
        g_v[i] = tagged(v, tag);
#endif
    }
}

} // namespace boltzmann_kernel

#ifdef USE_VERLET_DSFUN
typedef dsfloat float_type;
#else
typedef float float_type;
#endif

/**
 * device symbol wrappers
 */
// cuda::symbol<uint48>
//     boltzmann_wrapper<>::rand48::a = rng::gpu::rand48_kernel::a;
// cuda::symbol<uint48>
//     boltzmann_wrapper<>::rand48::c = rng::gpu::rand48_kernel::c;
// cuda::symbol<ushort3*>
//     boltzmann_wrapper<>::rand48::state = rng::gpu::rand48_kernel::g_state;

/**
 * device function wrappers
 */

cuda::function<void (float4*, uint, uint, float, float4*)>
    boltzmann_wrapper<3>::gaussian(boltzmann_kernel::gaussian<vector<float_type, 3> >);
cuda::function<void (float4*, uint, uint, float4 const*, dsfloat*)>
    boltzmann_wrapper<3>::shift_velocity(boltzmann_kernel::shift_velocity<vector<float_type, 3> >);
cuda::function<void (float4*, uint, uint, dsfloat const*, dsfloat)>
    boltzmann_wrapper<3>::scale_velocity(boltzmann_kernel::scale_velocity<vector<float_type, 3> >);

cuda::function<void (float4*, uint, uint, float, float2*)>
    boltzmann_wrapper<2>::gaussian(boltzmann_kernel::gaussian<vector<float_type, 2> >);
cuda::function<void (float4*, uint, uint, float2 const*, dsfloat*)>
    boltzmann_wrapper<2>::shift_velocity(boltzmann_kernel::shift_velocity<vector<float_type, 2> >);
cuda::function<void (float4*, uint, uint, dsfloat const*, dsfloat)>
    boltzmann_wrapper<2>::scale_velocity(boltzmann_kernel::scale_velocity<vector<float_type, 2> >);

}}} // namespace mdsim::gpu::velocity

} // namespace halmd
