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

#include <halmd/mdsim/gpu/base.cuh>
#include <halmd/mdsim/gpu/ljfluid_square.hpp>

using namespace halmd::gpu;

namespace halmd { namespace cu { namespace ljfluid
{

/**
 * MD simulation step
 */
template <typename vector_type,
          mixture_type mixture,
          potential_type potential,
          typename T>
__global__ void mdstep(float4 const* g_r, T* g_v, T* g_f, float* g_en, T* g_virial)
{
    enum { dimension = vector_type::static_size };

    extern __shared__ unsigned int s_tag[];
    vector_type* const s_r = reinterpret_cast<vector_type*>(&s_tag[TDIM]);

    // load particle associated with this thread
    unsigned int tag;
    vector_type r = detach_particle_tag(g_r[GTID], tag);
    vector_type v = g_v[GTID];
    // particle type in binary mixture
    int const a = (tag >= mpart[0]);

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    vector<float, (dimension - 1) * dimension / 2 + 1> virial = 0;
    // force sum
#ifdef USE_FORCE_DSFUN
    vector<dsfloat, dimension> f = 0;
#else
    vector_type f = 0;
#endif

    // iterate over all blocks
    for (unsigned int k = 0; k < gridDim.x; k++) {
        // load positions of particles within block
        __syncthreads();
        s_r[TID] = detach_particle_tag(g_r[k * blockDim.x + TID], s_tag[TID]);
        __syncthreads();

        // iterate over all particles within block
        for (unsigned int j = 0; j < blockDim.x; j++) {
            // skip placeholder particles
            if (k * blockDim.x + j >= npart)
                continue;
            // skip identical particle
            if (blockIdx.x == k && TID == j)
                continue;

            // particle type in binary mixture
            int const b = (s_tag[j] >= mpart[0]);
            // compute Lennard-Jones force with particle
            compute_force<mixture, potential>(r, s_r[j], f, en, virial, a + b);
        }
    }

    // second leapfrog step of integration of equations of motion
    leapfrog_full_step(v, static_cast<vector_type>(f));

    // store particle associated with this thread
    g_v[GTID] = v;
    g_f[GTID] = static_cast<vector_type>(f);
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

/**
 * sample trajectories
 */
template <int dimension, typename T>
__global__ void sample(float4 const* g_ir, T const* g_iR, T const* g_iv, T* g_or, T* g_ov)
{
    // permute particle phase space coordinates
    vector<float, dimension> const r = g_ir[GTID];
    vector<float, dimension> const R = g_iR[GTID];
    g_or[GTID] = r + box * R;
    g_ov[GTID] = g_iv[GTID];
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <int dimension, typename T>
__global__ void inteq(float4* g_r, T* g_R, T* g_v, T const* g_f)
{
    vector<float, dimension> r, dr, R, v, f;
    unsigned int tag;
    r = detach_particle_tag(g_r[GTID], tag);
    R = g_R[GTID];
    v = g_v[GTID];
    f = g_f[GTID];

    leapfrog_half_step(r, dr, R, v, f);

    g_r[GTID] = attach_particle_tag(r, tag);
    g_R[GTID] = R;
    g_v[GTID] = v;
}

}}} // namespace halmd::gpu::ljfluid

namespace halmd { namespace gpu
{

typedef ljfluid_base<ljfluid_impl_gpu_square> _Base;
typedef ljfluid<ljfluid_impl_gpu_square, 3> _3D;
typedef ljfluid<ljfluid_impl_gpu_square, 2> _2D;

/**
 * device function wrappers
 */
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<UNARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C0POT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<UNARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C2POT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<BINARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C0POT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<BINARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C2POT>);

cuda::function<void (float4 const*, float4 const*, float4 const*, float4*, float4*)>
    _3D::sample(cu::ljfluid::sample<3>);
cuda::function<void (float4*, float4*, float4*, float4 const*)>
    _3D::inteq(cu::ljfluid::inteq<3>);

cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<UNARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C0POT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<UNARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C2POT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<BINARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C0POT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<BINARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C2POT>);

cuda::function<void (float4 const*, float2 const*, float2 const*, float2*, float2*)>
    _2D::sample(cu::ljfluid::sample<2>);
cuda::function<void (float4*, float2*, float2*, float2 const*)>
    _2D::inteq(cu::ljfluid::inteq<2>);

}} // namespace halmd::gpu
