/* Maxwell-Boltzmann distribution at accurate temperature
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/int.hpp>
#include <ljgpu/algorithm/gpu/base.cuh>
#include <ljgpu/algorithm/gpu/reduce.cuh>
#include <ljgpu/math/gpu/dsfun.cuh>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>
#include <ljgpu/mdsim/gpu/boltzmann.hpp>
#define CU_NAMESPACE boltzmann
#include <ljgpu/rng/gpu/rand48.cuh>
using namespace boost;

namespace ljgpu { namespace cu { namespace boltzmann
{

enum { BLOCKS = ljgpu::gpu::boltzmann<>::BLOCKS };
enum { THREADS = ljgpu::gpu::boltzmann<>::THREADS };

/**
 * generate Maxwell-Boltzmann distributed velocities and reduce velocity
 */
template <typename vector_type, typename T>
__global__ void gaussian(T* g_v, uint npart, float temp, T* g_vcm)
{
    __shared__ vector_type s_vcm[THREADS];
    vector_type vcm = 0;

    for (uint i = GTID; i < npart; i += GTDIM) {
	T v;
	rand48::gaussian(v, temp);
	g_v[i] = v;
	vcm += v;
    }
    // reduced value for this thread
    s_vcm[TID] = vcm;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, sum_>(vcm, s_vcm);

    if (TID < 1) {
	// store block reduced value in global memory
	g_vcm[blockIdx.x] = vcm;
    }
}

/**
 * set center of mass velocity to zero and reduce squared velocity
 */
template <typename vector_type, typename T>
__global__ void shift_velocity(T* g_v, uint npart, T const* g_vcm, dfloat* g_vv)
{
    __shared__ vector_type s_vcm[BLOCKS];
    __shared__ dfloat s_vv[THREADS];
    vector_type vcm = 0;
    dfloat vv = 0;

    // compute mean center of mass velocity from block reduced values
    for (uint i = TID; i < BLOCKS; i += TDIM) {
	s_vcm[i] = g_vcm[i];
    }
    __syncthreads();
    for (uint i = 0; i < BLOCKS; ++i) {
	vcm += s_vcm[i];
    }
    vcm /= npart;

    for (uint i = GTID; i < npart; i += GTDIM) {
	vector_type v = g_v[i];
	v -= vcm;
	g_v[i] = v;
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
__global__ void scale_velocity(T* g_v, uint npart, dfloat const* g_vv, float temp)
{
    __shared__ dfloat s_vv[THREADS];
    dfloat vv = 0;

    // compute squared velocity sum from block reduced values
    for (uint i = TID; i < BLOCKS; i += TDIM) {
	s_vv[i] = g_vv[i];
    }
    __syncthreads();
    for (uint i = 0; i < BLOCKS; ++i) {
	vv += s_vv[i];
    }

    int dim = vector_type::static_size;
    float s = sqrtf(temp * dim * (npart / static_cast<float>(vv)));

    for (uint i = GTID; i < npart; i += GTDIM) {
	vector_type v = g_v[i];
	v *= s;
	g_v[i] = v;
    }
}

}}} // namespace ljgpu::cu::boltzmann

namespace ljgpu { namespace gpu
{

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
cuda::function<void (float4*, uint, float, float4*)>
    boltzmann<3>::gaussian(cu::boltzmann::gaussian<cu::vector<float, 3> >);
cuda::function<void (float4*, uint, float4 const*, dfloat*)>
    boltzmann<3>::shift_velocity(cu::boltzmann::shift_velocity<cu::vector<float, 3> >);
cuda::function<void (float4*, uint, dfloat const*, float)>
    boltzmann<3>::scale_velocity(cu::boltzmann::scale_velocity<cu::vector<float, 3> >);

cuda::function<void (float2*, uint, float, float2*)>
    boltzmann<2>::gaussian(cu::boltzmann::gaussian<cu::vector<float, 2> >);
cuda::function<void (float2*, uint, float2 const*, dfloat*)>
    boltzmann<2>::shift_velocity(cu::boltzmann::shift_velocity<cu::vector<float, 2> >);
cuda::function<void (float2*, uint, dfloat const*, float)>
    boltzmann<2>::scale_velocity(cu::boltzmann::scale_velocity<cu::vector<float, 2> >);

}} // namespace ljgpu::gpu
