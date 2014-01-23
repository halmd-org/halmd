/*
 * Copyright © 2011-2012  Michael Kopp
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
#include <halmd/mdsim/gpu/box_kernel.cuh> // reduce_periodic
#include <halmd/mdsim/gpu/mobilities/oseen_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh> // tagged/untagged
#include <halmd/numeric/mp/dsfloat.hpp> // sqrt
#include <halmd/utility/gpu/thread.cuh> // TID etc.

using namespace halmd::algorithm::gpu;
using namespace halmd::mdsim::gpu::particle_kernel;

//
// Compute velocities from forces via Oseen/Rotne-Prager tensor
//

namespace halmd {
namespace mdsim {
namespace gpu {
namespace mobilities {
namespace oseen_kernel {

/**
  * \brief Compute interaction mobility.
  *
  * @tparam order order of precision in (r/a) (1,2: oseen, >3: rotne-prager)
  * \note It is passed as template parameter, so that the compiler can decide
  * whether to implement oseen or rotne-prager part.
  */
template<
    int order
  , typename vector_type
  , typename vector_type_
>
__device__ void compute_mobility(
    vector_type const& r2
  , vector_type const& f2
  , vector_type const& r1
  , vector_type_& v1
  , vector_type const& box_length
  , float const radius
)
{
    vector_type dr = r1 - r2 ;

    // apply minimum image convention
    box_kernel::reduce_periodic(dr, box_length);

    float dr2 = inner_prod(dr, dr);
    float dr_norm = sqrtf(dr2);
    float b = radius / dr_norm;

    // to the actual oseen stuff
    if( order <= 2 ) { //oseen
        v1 += (f2 + (inner_prod(dr,f2) / dr2) * dr) * 0.75f * b;
    }
    else if (order <= 4) { // rotne prager
        if( dr_norm < 2*radius ) { // close branch
            v1 += ( 1 - (9.f / 32) * dr_norm / radius ) * f2 + ( (3.f / 32) * inner_prod(dr, f2) / (radius * dr_norm) ) * dr;
        }
        else { // default branch
            float b2 = b * b;
            v1 += ((0.75f + 0.5f * b2) * b) * f2 + ((0.75f - 1.5f * b2) * b * inner_prod(dr, f2) / dr2) * dr;
        }
    }
}


/**
  * \brief update velocities from positions using oseen tensor calculus
  *
  * Every thread computes velocity of one single (associated) particle;
  * thread GTID is responsible for which one (g_v[GTID]).
  *
  * The variable USE_OSEEN_DSFUN controls, whether the
  * computations internally are performed with single (if
  * not defined) or double (if defined) precision.  The velocity
  * is summed up from contributions from each particle --
  * USE_OSEEN_DSFUN controls, to which precision this sum is
  * evaluated.  The variable USE_VERLET_DSFUN on the other hand
  * controls, to what precision the velocities are stored.  Even
  * if it's not set (so velocities are stored in single
  * precision) it makes sense, to sum up the velocity in double
  * precision, and then store this sum in single precision.
  *
  * @param g_r positions in global device momory
  * @param g_f forces in global device momory
  * @param g_v velocities in global device momory -- will be updated in this function!
  * @param npart number of particles
  * @param radius hydrodynamic radius
  * @param self_mobility 1/(6*pi*eta*a) with eta being viscosity and a being radius
  * @tparam order order of precision in (r/a) (1,2: oseen, >3: rotne-prager)
  * @tparam vector_type float-vector type with appropriate dimension
  * @tparam vector_type_ dsfloat-vector type with appropriate dimension. If USE_OSEEN_DSFUN is set: dsfun. Else: float.
  * @tparam gpu_vector_type either float4 in 3D or float2 in 2D. Enables coalesced storage of forces.
  *
  * \note
  * In order to decrease the necessity for threads to read data
  * from the global memory, use shared memory of each block.  In
  * each step, the threads of one warp transfer data from global
  * to shared memory.  All threads of this block (which solely
  * have access to the shared memory) then read this data and
  * compute the velocities for `their' particle.  Then data from
  * the next particles are read -- again by threads of one warp.
  *
  * \note
  * In order not to use too much shared memory (especially for
  * older cards this is an important issue), only the threads of
  * one warp read from global memory.  As the card cannot handle
  * global read requests from different warps, this scheme is as
  * fast as if all warps would transfer data from global to
  * shared memory (as proposed in eg. CUDA C Programming Guide)
  * -- but it uses far less memory.
  *
  * \note
  * An alternative, and more simple, method is, to read all data
  * from global memory directly.  This is supposed to be slow
  * (cf. eg. CUDA C Programming Guide), but since the data for
  * is aligned so well in global memory, it's in fact fast (even
  * slightly faster than this algorihm).  Yet we have decided to
  * use the algorithm above, as it turned out to be quite stable
  * (even when data should not be aligned that well in global
  * memory it was fast) and we suppose, that it will adapt
  * better to new CUDA architectures.
  *
  */
template<
    int order
  , typename vector_type      // float-vector
  , typename vector_type_     // dsfun-vector
  , typename gpu_vector_type  // forces: gpu::vector<gpu_vector_type> (float4 in 3D, float2 in 2d)
>
__global__ void _compute_velocities(
    const float4* g_r
  , const gpu_vector_type* g_f
  , float4* g_v
  , const unsigned int npart
  , const vector_type box_length
  , const float radius
  , const float self_mobility
)
{
    // get information which thread this is and thus which particles are to be processed
    unsigned int const i = GTID;             // thread ID within grid
    unsigned int const threads_grid = GTDIM; // threads per grid

    // shared memory for this block
    //
    // CUDA only allows one pointer to shared memory.  Yet there is the
    // special construct
    //     extern __shared__ type name[];
    // which will create a (the) pointer to shared memory.  For this to work, a
    // default-shared-size must be passed to CUDA.  This is done via the third, optional,
    // parameter in cuda::configure(..).
    extern __shared__ char s_mem[];
    // position of other particles in shared memory
    float4* const s_r = reinterpret_cast<float4*>(s_mem);
    // forces of other particles in shared memory
    gpu_vector_type* const s_f = reinterpret_cast<gpu_vector_type*>(&s_r[WARP_SIZE]);

    // position of particle associated with this particular thread (single precision)
    //
    // Although for particles with i >= npart the following does not make sense
    // (as there are no positions to be fetched), it does not harm though.  The
    // particle module creates vectors big enough so that this operation will
    // not fail.  So for each thread connected to a ghost particle there will be
    // one superfluous access to global memory.  However if there was an
    // if-statement [if(i < npart)] before this, there would be one superflous
    // if statement for each single (real) particle.  So as there are
    // (normally) much more real than ghost particles, it makes sense to
    // simply apply these operations to the ghost ones, too...
    //
    // Similar situations in this file will be denoted by a `[*]'-symbol.
    vector_type r1 = g_r[i];

    // velocity of particle associated with this particular thread [*]
    vector_type_ v1;
    unsigned int this_tag;
#ifdef USE_OSEEN_DSFUN
#ifdef USE_VERLET_DSFUN
    tie(v1, this_tag) = untagged<vector_type_>(g_v[i], g_v[i + threads_grid]);
#else // ! USE_VERLET_DSFUN
    // vector_type_ is dsfun, but velocities are stored only in single precision.
    // Read zeros for all high precision bits.
    float4 _zero = {0, 0, 0, 0};
    tie(v1, this_tag) = untagged<vector_type_>(g_v[i], _zero);
#endif // end USE_VERLET_DSFUN
#else // ! USE_OSEEN_DSFUN
    tie(v1, this_tag) = untagged<vector_type_>(g_v[i]);
#endif // end USE_OSEEN_DSFUN
    // loop over every particle and consecutively add up velocity of this particle
    for (unsigned int offset = 0; offset < threads_grid; offset+=WARP_SIZE) {
        // transfer positions and forces from global to shared memory
        if (TID < WARP_SIZE) { // this is done by the first warp alone
            s_r[TID] = g_r[offset + TID];
            s_f[TID] = g_f[offset + TID];
        } // else { wait for this warp }

        // IMPORTANT: sync after reading. Otherwise a thread could request information not yet stored in shared memory
        __syncthreads();

        // loop over all data in shared memory
        for (unsigned int k = 0; k < WARP_SIZE; ++k ) {
            if (offset + k < npart) { //IMPORTANT: this must not be removed!
                // force on other particle
                vector_type f2 = s_f[k];

                if (i == offset+k) { // self mobility
                    v1 += f2;
                }
                else { // interaction
                    vector_type r2 = s_r[k];
                    compute_mobility<order>(r2, f2, r1, v1, box_length, radius);
                }
            }
        }
        __syncthreads(); // IMPORTANT: sync after computations
    }

    v1 *= self_mobility; // this has been factorized in previous computations

    // store final velocity for this particle [*]
#ifdef USE_OSEEN_DSFUN
#ifdef USE_VERLET_DSFUN
    tie(g_v[i], g_v[i + threads_grid]) = tagged(v1, this_tag);
#else // ! USE_VERLET_DSFUN
    // store velocity only in single precision
    // _zero is defined above and is overwritten...
    tie(g_v[i], _zero) = tagged(v1, this_tag);
#endif // end USE_VERLET_DSFUN
#else // ! USE_OSEEN_DSFUN
    g_v[i] = tagged(v1, this_tag);
#endif // end USE_OSEEN_DSFUN

}

} // namespace oseen_kernel

template <int dimension>
oseen_wrapper<dimension> const oseen_wrapper<dimension>::wrapper = {
#ifdef USE_OSEEN_DSFUN
    oseen_kernel::_compute_velocities<1, fixed_vector<float, dimension>, fixed_vector<dsfloat, dimension> > // _oseen
  , oseen_kernel::_compute_velocities<3, fixed_vector<float, dimension>, fixed_vector<dsfloat, dimension> > // _rotne
#else
    oseen_kernel::_compute_velocities<1, fixed_vector<float, dimension>, fixed_vector<float, dimension> > // _oseen
  , oseen_kernel::_compute_velocities<3, fixed_vector<float, dimension>, fixed_vector<float, dimension> > // _rotne
#endif
};

template class oseen_wrapper<3>;
template class oseen_wrapper<2>;

} // namespace mobilities
} // namespace gpu
} // namespace mdsim
} // namespace halmd
