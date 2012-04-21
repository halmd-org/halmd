/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/gpu/positions/lattice_kernel.hpp>
#include <halmd/mdsim/positions/lattice_primitive.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace positions {
namespace lattice_kernel {

template <typename vector_type, typename primitive_type>
__global__ void lattice(
    float4* g_r
  , unsigned int npart
  , float a
  , unsigned int skip
  , fixed_vector<float, vector_type::static_size> offset
  , fixed_vector<unsigned int, vector_type::static_size> ncell
)
{
    enum { dimension = vector_type::static_size };
    unsigned int const threads = GTDIM;

    for (uint i = GTID; i < npart; i += threads) {

        // load particle type
        vector_type r;
        unsigned int type;
#ifdef USE_VERLET_DSFUN
        tie(r, type) <<= tie(g_r[i], g_r[i + threads]);
#else
        tie(r, type) <<= g_r[i];
#endif

        // introduce a vacancy after every (skip - 1) particles
        uint nvacancies = (skip > 1) ? (i / (skip - 1)) : 0;

        // compute primitive lattice vector
        fixed_vector<float, dimension> e;
        primitive_type()(e, ncell, i + nvacancies);

        // scale with lattice constant and shift origin of lattice to offset
        r = e * a + offset; //< cast sum to dsfloat-based type

#ifdef USE_VERLET_DSFUN
        tie(g_r[i], g_r[i + threads]) <<= tie(r, type);
#else
        g_r[i] <<= tie(r, type);
#endif
    }
}

} // namespace lattice_kernel

template <int dimension>
lattice_wrapper<dimension> const lattice_wrapper<dimension>::kernel = {
#ifdef USE_VERLET_DSFUN
    lattice_kernel::lattice<fixed_vector<dsfloat, dimension>, fcc_lattice_primitive>
  , lattice_kernel::lattice<fixed_vector<dsfloat, dimension>, sc_lattice_primitive>
#else
    lattice_kernel::lattice<fixed_vector<float, dimension>, fcc_lattice_primitive>
  , lattice_kernel::lattice<fixed_vector<float, dimension>, sc_lattice_primitive>
#endif
};

template class lattice_wrapper<3>;
template class lattice_wrapper<2>;

} // namespace mdsim
} // namespace gpu
} // namespace positions
} // namespace halmd
