/*
 * Copyright © 2016       Manuel Dibak
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/gpu/orientations/uniform_kernel.hpp>
#include <halmd/mdsim/positions/lattice_primitive.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace orientations {
namespace uniform_kernel {

template <typename vector_type>
__global__ void uniform(
    float4* g_u
  , unsigned int npart
)
{
    enum { dimension = vector_type::static_size };
    unsigned int const threads = GTDIM;

    for (unsigned int i = GTID; i < npart; i += threads) {

        // load particle orientation
        vector_type u;
        unsigned int nothing;
#ifdef USE_VERLET_DSFUN
        tie(u, nothing) <<= tie(g_u[i], g_u[i + threads]);
#else
        tie(u, nothing) <<= g_u[i];
#endif
        u[0] = 1;
        u[1] = 0;
        u[2] = 0;

#ifdef USE_VERLET_DSFUN
        tie(g_u[i], g_u[i + threads]) <<= tie(u, nothing);
#else
        g_u[i] <<= tie(u, nothing);
#endif
    }
}

} // namespace uniform_kernel

uniform_wrapper const uniform_wrapper::kernel = {
#ifdef USE_VERLET_DSFUN
    uniform_kernel::uniform<fixed_vector<dsfloat, 3> >
#else
    uniform_kernel::uniform<fixed_vector<float, 3> >
#endif
};

//template class lattice_wrapper<close_packed_lattice<fixed_vector<float, 2>, fixed_vector<unsigned int, 2> > >;

} // namespace mdsim
} // namespace gpu
} // namespace orientations
} // namespace halmd
