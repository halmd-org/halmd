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

#ifndef HALMD_MDSIM_GPU_NEIGHBOR_WRAPPER_CUH
#define HALMD_MDSIM_GPU_NEIGHBOR_WRAPPER_CUH

#include <boost/mpl/if.hpp>

#include <cuda_wrapper.hpp>

namespace halmd { namespace mdsim { namespace gpu
{

template <size_t N>
struct neighbor_wrapper
{
    typedef typename boost::mpl::if_c<N == 3, float3, float2>::type vector_type;
    typedef typename boost::mpl::if_c<N == 3, uint3, uint2>::type cell_index;

    /** (cutoff lengths + neighbor list skin)² */
    static cuda::texture<float> rr_cut_skin;
    /** number of cells per dimension */
    static cuda::symbol<cell_index> ncell;
    /** neighbor list length */
    static cuda::symbol<unsigned int> neighbor_size;
    /** neighbor list stride */
    static cuda::symbol<unsigned int> neighbor_stride;
    /** number of particles in simulation box */
    static cuda::symbol<unsigned int> nbox;
    /** positions, tags */
    static cuda::texture<float4> r;
    /** cubic box edgle length */
    static cuda::symbol<vector_type> box_length;
    /** cell edge lengths */
    static cuda::symbol<vector_type> cell_length;
    /** assign particles to cells */
    static cuda::function<void (unsigned int*, unsigned int const*, unsigned int const*, unsigned int const*, unsigned int*)> assign_cells;
    /** compute global cell offsets in particle list */
    static cuda::function<void (unsigned int*, unsigned int*)> find_cell_offset;
    /** generate ascending index sequence */
    static cuda::function<void (unsigned int*)> gen_index;
    /** update neighbour lists */
    static cuda::function<void (unsigned int*, unsigned int*, unsigned int const*)> update_neighbours;
    /** compute cell indices for particle positions */
    static cuda::function<void (float4 const*, unsigned int*)> compute_cell;
};

}}} // namespace halmd::mdsim::gpu

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOR_WRAPPER_CUH */
