/*
 * Copyright © 2021      Jaslo Ziska
 * Copyright © 2008-2011 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/neighbours/from_particle_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {
namespace from_particle_kernel {

template <bool unroll_force_loop, typename vector_type>
__global__ void update(
    cudaTextureObject_t t_rr_cut_skin
  , float4 const* g_r1
  , unsigned int npart1
  , float4 const* g_r2
  , unsigned int npart2
  , unsigned int ntype1
  , unsigned int ntype2
  , vector_type box_length
  , unsigned int* g_neighbour
  , unsigned int size
  , unsigned int stride
  , int* g_overflow
)
{
    extern __shared__ unsigned int s_type[];
    vector_type* const s_r = reinterpret_cast<vector_type*>(&s_type[TDIM]);

    // load particle associated with this thread
    unsigned int index1 = GTID;
    unsigned int type1;
    vector_type r1;
    tie(r1, type1) <<= g_r1[index1];

    // number of particles in neighbour list
    unsigned int count = 0;

    // iterate over all blocks
    unsigned int nblock = (npart2 + blockDim.x - 1) / blockDim.x;
    for (unsigned int block = 0; block < nblock; ++block) {
        // load positions of particles within block
        __syncthreads();
        unsigned int index = block * blockDim.x + TID;
        if (index < npart2) {
            tie(s_r[TID], s_type[TID]) <<= g_r2[index];
        }
        __syncthreads();

        // skip placeholder particles
        if (index1 >= npart1) { continue; }

        // iterate over all particles within block
        for (unsigned int thread = 0; thread < blockDim.x; ++thread) {
            unsigned int index2 = block * blockDim.x + thread;
            // skip placeholder particles
            if (index2 >= npart2) { continue; }
            // skip identical particle
            if (g_r1 == g_r2 && index1 == index2) { continue; }

            vector_type r2 = s_r[thread];
            unsigned int type2 = s_type[thread];

            // particle distance vector
            vector_type dr = r1 - r2;
            // enforce periodic boundary conditions
            box_kernel::reduce_periodic(dr, box_length);
            // squared particle distance
            float rr = inner_prod(dr, dr);

            // enforce cutoff length with neighbour list skin
            float rr_cut_skin = tex1Dfetch<float>(t_rr_cut_skin, type1 * ntype2 + type2);

            if (rr > rr_cut_skin) { continue; }

            if (count < size) {
                // scattered write to neighbour list
                if (unroll_force_loop) {
                    g_neighbour[index1 * size + count] = index2;
                } else {
                    g_neighbour[count * stride + index1] = index2;
                }
                // increment neighbour list particle count
                count++;
            }
            else {
                atomicAdd(g_overflow, 1);
            }
        }
    }
}

template<typename vector_type>
int block_size_to_smem_size(int block_size) {
    return block_size * (sizeof(unsigned int) + sizeof(vector_type));
}

} // namespace from_particle_kernel

template <int dimension>
from_particle_wrapper<dimension> from_particle_wrapper<dimension>::kernel = {
    from_particle_kernel::update<true, fixed_vector<float, dimension>>
  , from_particle_kernel::update<false, fixed_vector<float, dimension>>
};

template class from_particle_wrapper<3>;
template class from_particle_wrapper<2>;

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd
