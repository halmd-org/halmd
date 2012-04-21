/*
 * Copyright © 2011-2012  Michael Kopp and Felix Höfling
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

#include <halmd/mdsim/gpu/integrators/euler_kernel.cuh>
#include <halmd/mdsim/gpu/integrators/euler_kernel.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {
namespace euler_kernel {

/**
 * Euler integration
 *
 * @param g_r           positions
 * @param g_image       number of times the particle exceeded the box margin
 * @param g_v           velocities
 * @param g_f           forces
 * @param timestep      integration timestep
 * @param box_length    edge lengths of cuboid box
 */
template <
    typename vector_type //< precision of positions and velocities
  , typename vector_type_ //< default precision (box size, image vector)
  , typename gpu_vector_type
>
__global__ void _integrate(
    float4* g_r
  , gpu_vector_type* g_image
  , float4* g_v
  , float timestep
  , vector_type_ box_length
)
{
    // get information which thread this is and thus which particles are to
    // be processed
    unsigned int const i = GTID;
    unsigned int const threads = GTDIM;
    unsigned int type, tag;
    // local copy of position and velocity
    vector_type r, v;
#ifdef USE_VERLET_DSFUN
    tie(r, type) <<= tie(g_r[i], g_r[i + threads]);
    tie(v, tag) <<= tie(g_v[i], g_v[i + threads]);
#else
    tie(r, type) <<= g_r[i];
    tie(v, tag) <<= g_v[i];
#endif
    // local copy of image
    vector_type_ image = g_image[i];

    // run actual integration routine in .cuh file
    integrate(r, image, v, timestep, box_length);

#ifdef USE_VERLET_DSFUN
    tie(g_r[i], g_r[i + threads]) <<= tie(r, type);
#else
    g_r[i] <<= tie(r, type);
#endif
    g_image[i] = image;
}

} // namespace euler_kernel

template <int dimension>
euler_wrapper<dimension> const euler_wrapper<dimension>::kernel = {
#ifdef USE_VERLET_DSFUN
    euler_kernel::_integrate<fixed_vector<dsfloat, dimension> >
#else
    euler_kernel::_integrate<fixed_vector<float, dimension> >
#endif
};

// explicit instantiation
template class euler_wrapper<3>;
template class euler_wrapper<2>;

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
