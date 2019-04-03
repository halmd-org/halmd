/*
 * Copyright © 2014-2015 Nicolas Höft
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

#include <cub/iterator/counting_input_iterator.cuh>

#include <halmd/algorithm/gpu/copy_if_kernel.cuh>
#include <halmd/mdsim/geometries/cuboid.hpp>
#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/region_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace region_kernel {

template <typename vector_type, typename geometry_type>
__global__ void compute_mask(
    float4 const* g_r
  , unsigned int nparticle
  , unsigned int* g_mask
  , geometry_type const geometry
  , geometry_selection selection
  , vector_type box_length
)
{
    enum { dimension = vector_type::static_size };
    unsigned int const i = GTID;
    if(i >= nparticle)
        return;

    vector_type r;
    unsigned int type;
    tie(r, type) <<= g_r[i];

    // enforce periodic boundary conditions
    box_kernel::reduce_periodic(r, box_length);
    bool in_geometry = geometry(r);
    if(selection == excluded)
        in_geometry = !in_geometry;
    // 1 means the particle in in the selector, 0 means outside
    g_mask[i] = in_geometry ? 1 : 0;
}

template<typename geometry_type>
struct geometry_predicate
{
    typedef typename geometry_type::vector_type vector_type;

    geometry_predicate(float4 const* position, geometry_type const& geometry, geometry_selection sel)
      : position_(position)
      , geometry_(geometry)
      , selection_(sel)
    {};

    HALMD_GPU_ENABLED bool operator()(int const& i) const
    {
        vector_type r;
        unsigned int type;
        tie(r, type) <<= position_[i];
        bool in_geometry = geometry_(r);
        if (selection_ == excluded)
            return !in_geometry;
        return in_geometry;
    }

private:
    float4 const* position_; // position array
    geometry_type const geometry_;
    geometry_selection selection_;

};

template <typename geometry_type>
unsigned int copy_selection(float4 const* g_r, unsigned int nparticle, unsigned int* g_output, geometry_type const geometry, geometry_selection selection)
{
    typedef geometry_predicate<geometry_type> predicate_type;
    predicate_type const predicate(g_r, geometry, selection);

    // iterate over the particle indices, not the positions itself
    cub::CountingInputIterator<int> index(0);
    unsigned int output_size =
    halmd::algorithm::gpu::copy_if_kernel::copy_if(
        index
      , nparticle
      , predicate
      , g_output
    );

    return output_size;
}

} // namespace region_kernel

template<int dimension, typename geometry_type>
region_wrapper<dimension, geometry_type> const
region_wrapper<dimension, geometry_type>::kernel = {
    region_kernel::compute_mask
  , region_kernel::copy_selection<geometry_type>
};

template class region_wrapper<3, halmd::mdsim::geometries::cuboid<3, float> >;
template class region_wrapper<2, halmd::mdsim::geometries::cuboid<2, float> >;

} // namespace gpu
} // namespace mdsim

} // namespace halmd
