/* Lennard-Jones fluid kernel
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
#include <halmd/math/gpu/vector2d.cuh>
#include <halmd/math/gpu/vector3d.cuh>
#include <halmd/mdsim/gpu/lattice.hpp>

namespace halmd { namespace cu { namespace lattice
{

/**
 * place particles on a face centered cubic lattice (fcc)
 */
__device__ void fcc(vector<float, 3>& r, uint n)
{
    // compose primitive vectors from 1-dimensional index
    r.x = ((GTID >> 2) % n) + ((GTID ^ (GTID >> 1)) & 1) / 2.f;
    r.y = ((GTID >> 2) / n % n) + (GTID & 1) / 2.f;
    r.z = ((GTID >> 2) / n / n) + (GTID & 2) / 4.f;
}

__device__ void fcc(vector<float, 2>& r, uint n)
{
    r.x = ((GTID >> 1) % n) + (GTID & 1) / 2.f;
    r.y = ((GTID >> 1) / n) + (GTID & 1) / 2.f;
}

/**
 * place particles on a simple cubic lattice (sc)
 */
__device__ void sc(vector<float, 3>& r, uint n)
{
    r.x = (GTID % n) + 0.5f;
    r.y = (GTID / n % n) + 0.5f;
    r.z = (GTID / n / n) + 0.5f;
}

__device__ void sc(vector<float, 2>& r, uint n)
{
    r.x = (GTID % n) + 0.5f;
    r.y = (GTID / n) + 0.5f;
}

template <int dimension, void (*primitive)(vector<float, dimension>&, uint)>
__global__ void lattice(float4* g_r, uint n, float box)
{
    vector<float, dimension> r;
    primitive(r, n);
    g_r[GTID] = r * (box / n);
}

}}} // namespace halmd::cu::lattice

namespace halmd { namespace gpu
{

typedef lattice<3> __3D;
typedef lattice<2> __2D;

/**
 * device function wrappers
 */
cuda::function<void (float4*, uint, float)>
    __3D::fcc(cu::lattice::lattice<3, cu::lattice::fcc>);
cuda::function<void (float4*, uint, float)>
    __3D::sc(cu::lattice::lattice<3, cu::lattice::sc>);
cuda::function<void (float4*, uint, float)>
    __2D::fcc(cu::lattice::lattice<2, cu::lattice::fcc>);
cuda::function<void (float4*, uint, float)>
    __2D::sc(cu::lattice::lattice<2, cu::lattice::sc>);

}} // namespace halmd::gpu
