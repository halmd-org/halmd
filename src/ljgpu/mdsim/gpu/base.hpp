/* Lennard-Jones fluid kernel
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

#ifndef LJGPU_MDSIM_GPU_LJFLUID_BASE_HPP
#define LJGPU_MDSIM_GPU_LJFLUID_BASE_HPP

#include <cuda_wrapper.hpp>
#include <ljgpu/mdsim/impl.hpp>
#include <ljgpu/mdsim/variant.hpp>
#include <ljgpu/rng/gpu/uint48.cuh>

namespace ljgpu { namespace gpu
{

enum { VIRTUAL_PARTICLE = -1 };

/**
 * extract particle number from particle tag
 */
__device__ __host__ inline int particle_id(int tag)
{
    return (tag & 0x3FFFFFFF);
}

/*
 * extract particle type from particle tag
 */
__device__ __host__ inline int particle_type(int tag)
{
    return ((tag >> 30) & 0x1);
}

/**
 * compute particle tag from number and type
 */
__device__ __host__ inline int particle_tag(int id, int type = 0)
{
    return ((id & 0x3FFFFFFF) | ((type & 0x1) << 30));
}


template <template <int> class ljfluid_impl>
struct ljfluid_base;

template <>
struct ljfluid_base<ljfluid_impl_gpu_base>
{
    static cuda::symbol<uint> npart;
    static cuda::symbol<float> box;
    static cuda::symbol<float> timestep;
    static cuda::symbol<float[]> r_cut;
    static cuda::symbol<float[]> rr_cut;
    static cuda::symbol<float> en_cut;
    /** binary mixture */
    static cuda::symbol<uint[]> mpart;
    static cuda::symbol<float[]> epsilon;
    static cuda::symbol<float[]> sigma2;
    /** C² potential */
    static cuda::symbol<float> rri_smooth;
    /** NVT ensemble */
    static cuda::symbol<float> thermostat_nu;
    static cuda::symbol<float> thermostat_temp;

    struct rand48
    {
	static cuda::symbol<uint48> a;
	static cuda::symbol<uint48> c;
	static cuda::symbol<ushort3*> state;
    };

    static cuda::function<void (float3*, const float2)> sample_smooth_function;
};

template <typename ljfluid_impl>
struct ljfluid;

template <>
struct ljfluid<ljgpu::ljfluid_impl_gpu_base<3> >
: public ljfluid_base<ljfluid_impl_gpu_base>
{
    static cuda::function<void (float4*, float4*, float4*, float4 const*)> inteq;
    static cuda::function<void (float4*, float)> boltzmann;
    static cuda::function<void (float4*, int*)> init_tags;
    static cuda::function<void (float4*, int*)> init_types;
};

template <>
struct ljfluid<ljgpu::ljfluid_impl_gpu_base<2> >
: public ljfluid_base<ljfluid_impl_gpu_base>
{
    static cuda::function<void (float4*, float2*, float2*, float2 const*)> inteq;
    static cuda::function<void (float2*, float)> boltzmann;
    static cuda::function<void (float4*, int*)> init_tags;
    static cuda::function<void (float4*, int*)> init_types;
};

}} // namespace ljgpu::gpu

#endif /* ! LJGPU_MDSIM_GPU_LJFLUID_BASE_HPP */
