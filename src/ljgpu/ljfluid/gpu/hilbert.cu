/* Hilbert spacing-filling curve kernel
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

#include <float.h>
#include <ljgpu/algorithm/gpu/base.cuh>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>
#include <ljgpu/ljfluid/gpu/hilbert.hpp>
using namespace ljgpu::gpu::hilbert;

namespace ljgpu { namespace cu { namespace hilbert
{

/** periodic box length */
__constant__ float box;
/** Hilbert space-filling curve recursion depth */
__constant__ unsigned int depth;

/**
 * swap Hilbert spacing-filling curve vertices
 */
__device__ void vertex_swap(uint& v, uint& a, uint& b, uint const& mask)
{
    // swap bits comprising Hilbert codes in vertex-to-code lookup table
    const uint va = ((v >> a) & mask);
    const uint vb = ((v >> b) & mask);
    v = v ^ (va << a) ^ (vb << b) ^ (va << b) ^ (vb << a);
    // update code-to-vertex lookup table
    swap(a, b);
}

/**
 * map 3-dimensional point to 1-dimensional point on Hilbert space curve
 */
__global__ void hilbert_curve(float4 const* g_r, unsigned int* g_sfc)
{
    //
    // Jun Wang & Jie Shan, Space-Filling Curve Based Point Clouds Index,
    // GeoComputation, 2005
    //

    // Hilbert code for particle
    unsigned int hcode = 0;

    //
    // We need to avoid ambiguities during the assignment of a particle
    // to a subcell, i.e. the particle position should never lie on an
    // edge or corner of multiple subcells, or the algorithm will have
    // trouble converging to a definite Hilbert curve.
    //
    // Therefore, we use a simple cubic lattice of predefined dimensions
    // according to the number of cells at the deepest recursion level,
    // and round the particle position to the nearest center of a cell.
    //

    // Hilbert cells per dimension at deepest recursion level
    const uint n = 1UL << depth;
    // fractional index of particle's Hilbert cell in [0, n)
    const float3 cell = (__saturatef(unpack(g_r[GTID]) / box) * (1.f - FLT_EPSILON)) * n;

    // round particle position to center of cell in unit coordinates
    float3 r = (floorf(cell) + make_float3(0.5f, 0.5f, 0.5f)) / n;
    // use symmetric coordinates
    r -= make_float3(0.5f, 0.5f, 0.5f);

    // Hilbert code-to-vertex lookup table
    uint a = 21;
    uint b = 18;
    uint c = 12;
    uint d = 15;
    uint e = 3;
    uint f = 0;
    uint g = 6;
    uint h = 9;
    // Hilbert vertex-to-code lookup table
    uint vc = 1U << b ^ 2U << c ^ 3U << d ^ 4U << e ^ 5U << f ^ 6U << g ^ 7U << h;

#define MASK ((1 << 3) - 1)

    // 32-bit integer for 3D Hilbert code allows a maximum of 10 levels
    for (unsigned int i = 0; i < depth; ++i) {
	// determine Hilbert vertex closest to particle
	const uint x = __signbitf(r.x) & 1;
	const uint y = __signbitf(r.y) & 1;
	const uint z = __signbitf(r.z) & 1;
	// lookup Hilbert code
	const uint v = (vc >> (3 * (x + (y << 1) + (z << 2))) & MASK);

	// scale particle coordinates to subcell
	r = 2 * r - make_float3(0.5f - x, 0.5f - y, 0.5f - z);
	// apply permutation rule according to Hilbert code
	if (v == 0) {
	    vertex_swap(vc, b, h, MASK);
	    vertex_swap(vc, c, e, MASK);
	}
	else if (v == 1 || v == 2) {
	    vertex_swap(vc, c, g, MASK);
	    vertex_swap(vc, d, h, MASK);
	}
	else if (v == 3 || v == 4) {
	    vertex_swap(vc, a, c, MASK);
#ifdef USE_ALTERNATIVE_HILBERT_3D
	    vertex_swap(vc, b, d, MASK);
	    vertex_swap(vc, e, g, MASK);
#endif
	    vertex_swap(vc, f, h, MASK);
	}
	else if (v == 5 || v == 6) {
	    vertex_swap(vc, a, e, MASK);
	    vertex_swap(vc, b, f, MASK);
	}
	else if (v == 7) {
	    vertex_swap(vc, a, g, MASK);
	    vertex_swap(vc, d, f, MASK);
	}

	// add vertex code to partial Hilbert code
	hcode = (hcode << 3) + v;
    }

#undef MASK

    // store Hilbert code for particle
    g_sfc[GTID] = hcode;
}

__global__ void hilbert_curve(float2 const* g_r, unsigned int* g_sfc)
{
    // Hilbert code for particle
    unsigned int hcode = 0;
    // Hilbert cells per dimension at deepest recursion level
    const uint n = 1UL << depth;
    // fractional index of particle's Hilbert cell in [0, n)
    const float2 cell = (__saturatef(unpack(g_r[GTID]) / box) * (1.f - FLT_EPSILON)) * n;

    // round particle position to center of cell in unit coordinates
    float2 r = (floorf(cell) + make_float2(0.5f, 0.5f)) / n;
    // use symmetric coordinates
    r -= make_float2(0.5f, 0.5f);

    // Hilbert code-to-vertex lookup table
    uint a = 6;
    uint b = 4;
    uint c = 0;
    uint d = 2;
    // Hilbert vertex-to-code lookup table
    uint vc = 1U << b ^ 2U << c ^ 3U << d;

#define MASK ((1 << 2) - 1)

    // 32-bit integer for 2D Hilbert code allows a maximum of 16 levels
    for (unsigned int i = 0; i < depth; ++i) {
	// determine Hilbert vertex closest to particle
	const uint x = __signbitf(r.x) & 1;
	const uint y = __signbitf(r.y) & 1;
	// lookup Hilbert code
	const uint v = (vc >> (2 * (x + (y << 1))) & MASK);

	// scale particle coordinates to subcell
	r = 2 * r - make_float2(0.5f - x, 0.5f - y);
	// apply permutation rule according to Hilbert code
	if (v == 0) {
	    vertex_swap(vc, b, d, MASK);
	}
	else if (v == 3) {
	    vertex_swap(vc, a, c, MASK);
	}

	// add vertex code to partial Hilbert code
	hcode = (hcode << 2) + v;
    }

#undef MASK

    // store Hilbert code for particle
    g_sfc[GTID] = hcode;
}

}}} // namespace ljgpu::cu::hilbert

namespace ljgpu { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<void (float4 const*, unsigned int*),
	       void (float2 const*, unsigned int*)>
    hilbert::curve(cu::hilbert::hilbert_curve, cu::hilbert::hilbert_curve);

/**
 * device constant wrappers
 */
cuda::symbol<float> hilbert::box(cu::hilbert::box);
cuda::symbol<unsigned int> hilbert::depth(cu::hilbert::depth);

}} // namespace ljgpu::gpu
