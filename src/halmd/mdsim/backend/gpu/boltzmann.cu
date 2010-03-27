/* Maxwell-Boltzmann distribution at accurate temperature
 *
 * Copyright © 2008-2009  Peter Colberg
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

#include <halmd/algorithm/gpu/base.cuh>
#include <halmd/algorithm/gpu/reduce.cuh>
#include <halmd/math/gpu/dsvector.cuh>
#include <halmd/math/gpu/vector2d.cuh>
#include <halmd/math/gpu/vector3d.cuh>
#include <halmd/mdsim/backend/gpu/boltzmann.hpp>
#define CU_NAMESPACE boltzmann
#include <halmd/rng/gpu/rand48.cuh>

namespace halmd { namespace cu { namespace boltzmann
{

enum { BLOCKS = halmd::gpu::boltzmann<>::BLOCKS };
enum { THREADS = halmd::gpu::boltzmann<>::THREADS };

/**
 * generate Maxwell-Boltzmann distributed velocities and reduce velocity
 */
template <typename vector_type, typename T>
__global__ void gaussian(T* g_v, uint npart, uint nplace, float temp, T* g_vcm)
{
    enum { dimension = vector_type::static_size };
    __shared__ vector_type s_vcm[THREADS];
    vector_type vcm = 0;

    // read random number generator state from global device memory
    rand48::state_type state = rand48::g_state[GTID];

    for (uint i = GTID; i < npart; i += GTDIM) {
        T v;
        rand48::gaussian(v, temp, state);
        g_v[i] = v;
#ifdef USE_VERLET_DSFUN
        g_v[i + nplace] = vector<float, dimension>(0);
#endif
        vcm += vector<float, dimension>(v);
    }
    // store random number generator state in global device memory
    rand48::g_state[GTID] = state;

    // reduced value for this thread
    s_vcm[TID] = vcm;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, sum_>(vcm, s_vcm);

    if (TID < 1) {
        // store block reduced value in global memory
        g_vcm[blockIdx.x] = static_cast<vector<float, dimension> >(vcm);
#ifdef USE_VERLET_DSFUN
        g_vcm[blockIdx.x + BDIM] = dsfloat2lo(vcm);
#endif
    }
}

/**
 * set center of mass velocity to zero and reduce squared velocity
 */
template <typename vector_type, typename T>
__global__ void shift_velocity(T* g_v, uint npart, uint nplace, T const* g_vcm, dsfloat* g_vv)
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
        vector_type v(g_v[i], g_v[i + nplace]);
#else
        vector_type v = g_v[i];
#endif
        v -= vcm;
        g_v[i] = static_cast<vector<float, dimension> >(v);
#ifdef USE_VERLET_DSFUN
        g_v[i + nplace] = dsfloat2lo(v);
#endif
        vv += v * v;
    }
    // reduced value for this thread
    s_vv[TID] = vv;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, sum_>(vv, s_vv);

    if (TID < 1) {
        // store block reduced value in global memory
        g_vv[blockIdx.x] = vv;
    }
}

/**
 * rescale velocities to accurate temperature
 */
template <typename vector_type, typename T>
__global__ void scale_velocity(T* g_v, uint npart, uint nplace, dsfloat const* g_vv, dsfloat temp)
{
    enum { dimension = vector_type::static_size };
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

    int dim = vector_type::static_size;
    dsfloat coeff = sqrt(temp * static_cast<dsfloat>(dim) * (static_cast<dsfloat>(npart) / vv));

    for (uint i = GTID; i < npart; i += GTDIM) {
#ifdef USE_VERLET_DSFUN
        vector_type v(g_v[i], g_v[i + nplace]);
#else
        vector_type v = g_v[i];
#endif
        v *= coeff;
        g_v[i] = static_cast<vector<float, dimension> >(v);
#ifdef USE_VERLET_DSFUN
        g_v[i + nplace] = dsfloat2lo(v);
#endif
    }
}

}}} // namespace halmd::cu::boltzmann

namespace halmd { namespace gpu
{

#ifdef USE_VERLET_DSFUN
typedef dsfloat float_type;
#else
typedef float float_type;
#endif

/**
 * device symbol wrappers
 */
cuda::symbol<uint48>
    boltzmann<>::rand48::a(cu::boltzmann::rand48::a);
cuda::symbol<uint48>
    boltzmann<>::rand48::c(cu::boltzmann::rand48::c);
cuda::symbol<ushort3*>
    boltzmann<>::rand48::state(cu::boltzmann::rand48::g_state);

/**
 * device function wrappers
 */

cuda::function<void (float4*, uint, uint, float, float4*)>
    boltzmann<3>::gaussian(cu::boltzmann::gaussian<cu::vector<float_type, 3> >);
cuda::function<void (float4*, uint, uint, float4 const*, dsfloat*)>
    boltzmann<3>::shift_velocity(cu::boltzmann::shift_velocity<cu::vector<float_type, 3> >);
cuda::function<void (float4*, uint, uint, dsfloat const*, dsfloat)>
    boltzmann<3>::scale_velocity(cu::boltzmann::scale_velocity<cu::vector<float_type, 3> >);

cuda::function<void (float2*, uint, uint, float, float2*)>
    boltzmann<2>::gaussian(cu::boltzmann::gaussian<cu::vector<float_type, 2> >);
cuda::function<void (float2*, uint, uint, float2 const*, dsfloat*)>
    boltzmann<2>::shift_velocity(cu::boltzmann::shift_velocity<cu::vector<float_type, 2> >);
cuda::function<void (float2*, uint, uint, dsfloat const*, dsfloat)>
    boltzmann<2>::scale_velocity(cu::boltzmann::scale_velocity<cu::vector<float_type, 2> >);

}} // namespace halmd::gpu
