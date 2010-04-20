/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <boost/mpl/if.hpp>
#include <float.h>

#include <halmd/algorithm/gpu/base.cuh>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>

using namespace boost::mpl;
using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::numeric::gpu::blas;

namespace halmd { namespace mdsim { namespace gpu
{

namespace hilbert_kernel
{

template <size_t N>
struct dim_
{
    /** positions, tags */
    static texture<float4, 1, cudaReadModeElementType> r;
    /** cubic box edgle length */
    static __constant__ typename if_c<N == 3, float3, float2>::type box_length;
};

// explicit instantiation
template class dim_<3>;
template class dim_<2>;

/** Hilbert space-filling curve recursion depth */
__constant__ unsigned int depth_;

/**
 * swap Hilbert spacing-filling curve vertices
 */
__device__ void swap(uint& v, uint& a, uint& b, uint const& mask)
{
    // swap bits comprising Hilbert codes in vertex-to-code lookup table
    uint const va = ((v >> a) & mask);
    uint const vb = ((v >> b) & mask);
    v = v ^ (va << a) ^ (vb << b) ^ (va << b) ^ (vb << a);
    // update code-to-vertex lookup table
    algorithm::gpu::swap(a, b);
}

/**
 * map 3-dimensional point to 1-dimensional point on Hilbert space curve
 */
__device__ unsigned int _map(vector<float, 3> r)
{
    //
    // Jun Wang & Jie Shan, Space-Filling Curve Based Point Clouds Index,
    // GeoComputation, 2005
    //

    // Hilbert code for particle
    unsigned int hcode = 0;
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
    for (unsigned int i = 0; i < depth_; ++i) {
        // determine Hilbert vertex closest to particle
        vector<unsigned int, 3> x;
        x[0] = __signbitf(r[0]) & 1;
        x[1] = __signbitf(r[1]) & 1;
        x[2] = __signbitf(r[2]) & 1;
        // lookup Hilbert code
        const uint v = (vc >> (3 * (x[0] + (x[1] << 1) + (x[2] << 2))) & MASK);

        // scale particle coordinates to subcell
        r = 2 * r - (vector<float, 3>(0.5f) - vector<float, 3>(x));
        // apply permutation rule according to Hilbert code
        if (v == 0) {
            swap(vc, b, h, MASK);
            swap(vc, c, e, MASK);
        }
        else if (v == 1 || v == 2) {
            swap(vc, c, g, MASK);
            swap(vc, d, h, MASK);
        }
        else if (v == 3 || v == 4) {
            swap(vc, a, c, MASK);
#ifdef USE_HILBERT_ALT_3D
            swap(vc, b, d, MASK);
            swap(vc, e, g, MASK);
#endif
            swap(vc, f, h, MASK);
        }
        else if (v == 5 || v == 6) {
            swap(vc, a, e, MASK);
            swap(vc, b, f, MASK);
        }
        else if (v == 7) {
            swap(vc, a, g, MASK);
            swap(vc, d, f, MASK);
        }

        // add vertex code to partial Hilbert code
        hcode = (hcode << 3) + v;
    }
#undef MASK
    return hcode;
}

/**
 * map 2-dimensional point to 1-dimensional point on Hilbert space curve
 */
__device__ unsigned int _map(vector<float, 2> r)
{
    // Hilbert code for particle
    unsigned int hcode = 0;
    // Hilbert code-to-vertex lookup table
    uint a = 6;
    uint b = 4;
    uint c = 0;
    uint d = 2;
    // Hilbert vertex-to-code lookup table
    uint vc = 1U << b ^ 2U << c ^ 3U << d;

#define MASK ((1 << 2) - 1)

    // 32-bit integer for 2D Hilbert code allows a maximum of 16 levels
    for (unsigned int i = 0; i < depth_; ++i) {
        // determine Hilbert vertex closest to particle
        vector<unsigned int, 2> x;
        x[0] = __signbitf(r[0]) & 1;
        x[1] = __signbitf(r[1]) & 1;
        // lookup Hilbert code
        const uint v = (vc >> (2 * (x[0] + (x[1] << 1))) & MASK);

        // scale particle coordinates to subcell
        r = 2 * r - (vector<float, 2>(0.5f) - vector<float, 2>(x));
        // apply permutation rule according to Hilbert code
        if (v == 0) {
            swap(vc, b, d, MASK);
        }
        else if (v == 3) {
            swap(vc, a, c, MASK);
        }

        // add vertex code to partial Hilbert code
        hcode = (hcode << 2) + v;
    }
#undef MASK
    return hcode;
}

/**
 * generate Hilbert space-filling curve
 */
template <typename vector_type>
__global__ void map(float4 const* g_r, unsigned int* g_sfc)
{
    enum { dimension = vector_type::static_size };

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

    vector_type r = untagged<vector_type>(g_r[GTID]);
    vector_type L = dim_<dimension>::box_length;
    // Hilbert cells per dimension at deepest recursion level
    uint const n = 1UL << depth_;
    // fractional index of particle's Hilbert cell in [0, n)
    r = n * (__saturate(element_div(r, L)) * (1.f - FLT_EPSILON));

    // round particle position to center of cell in unit coordinates
    r = (floor(r) + vector_type(0.5f)) / n;
    // use symmetric coordinates
    r -= vector_type(0.5f);

    // compute Hilbert code for particle
    g_sfc[GTID] = _map(r);
}

} // namespace hilbert_kernel

}}} //namespace halmd::mdsim::gpu
