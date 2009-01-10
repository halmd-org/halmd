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

#include <ljgpu/algorithm/gpu/base.cuh>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>
#include <ljgpu/ljfluid/gpu/lattice.hpp>
using namespace ljgpu::gpu::lattice;

namespace ljgpu { namespace gpu
{

/**
 * place particles on a face centered cubic lattice (fcc)
 */
__global__ void fcc(float4* g_r, uint n, float box)
{
    float3 r;
    // compose primitive vectors from 1-dimensional index
    r.x = ((GTID >> 2) % n) + ((GTID ^ (GTID >> 1)) & 1) / 2.f;
    r.y = ((GTID >> 2) / n % n) + (GTID & 1) / 2.f;
    r.z = ((GTID >> 2) / n / n) + (GTID & 2) / 4.f;
    g_r[GTID] = pack(r * (box / n));
}

__global__ void fcc(float2* g_r, uint n, float box)
{
    float2 r;
    r.x = ((GTID >> 1) % n) + (GTID & 1) / 2.f;
    r.y = ((GTID >> 1) / n) + (GTID & 1) / 2.f;
    g_r[GTID] = pack(r * (box / n));
}

/**
 * place particles on a simple cubic lattice (sc)
 */
__global__ void sc(float4* g_r, uint n, float box)
{
    float3 r;
    r.x = (GTID % n) + 0.5f;
    r.y = (GTID / n % n) + 0.5f;
    r.z = (GTID / n / n) + 0.5f;
    g_r[GTID] = pack(r * (box / n));
}

__global__ void sc(float2* g_r, uint n, float box)
{
    float2 r;
    r.x = (GTID % n) + 0.5f;
    r.y = (GTID / n) + 0.5f;
    g_r[GTID] = pack(r * (box / n));
}

/**
 * device function wrappers
 */
cuda::function<void (float2*, uint, float), void (float4*, uint, float)>
    lattice::fcc(gpu::lattice_fcc, gpu::lattice_fcc);
cuda::function<void (float2*, uint, float), void (float4*, uint, float)>
    lattice::sc(gpu::lattice_sc, gpu::lattice_sc);

}} // namespace ljgpu::gpu
