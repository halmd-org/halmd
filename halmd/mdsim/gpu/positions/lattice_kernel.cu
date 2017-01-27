/*
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

template <typename ptr_type, typename vector_type, typename lattice_type>
__global__ void lattice(
    ptr_type g_r
  , unsigned int npart
  , float a
  , unsigned int skip
  , typename lattice_type::result_type offset
  , typename lattice_type::shape_type ncell
)
{
    enum { dimension = vector_type::static_size };
    unsigned int const threads = GTDIM;

    lattice_type const lattice(ncell);

    for (unsigned int i = GTID; i < npart; i += threads) {

        // load particle type
        vector_type r;
        unsigned int type;
        tie(r, type) <<= g_r[i];

        // introduce a vacancy after every (skip - 1) particles
        uint nvacancies = (skip > 1) ? (i / (skip - 1)) : 0;

        // compute primitive lattice vector
        fixed_vector<float, dimension> e = lattice(i + nvacancies);

        // scale with lattice constant and shift origin of lattice to offset
        r = e * a + offset; //< cast sum to dsfloat-based type

        g_r[i] <<= tie(r, type);
    }
}

} // namespace lattice_kernel

template <typename float_type, typename lattice_type>
lattice_wrapper<float_type, lattice_type> const lattice_wrapper<float_type, lattice_type>::kernel = {
    lattice_kernel::lattice<ptr_type, fixed_vector<float_type, dimension>, lattice_type>
};

template class lattice_wrapper<float, close_packed_lattice<fixed_vector<float, 3>, fixed_vector<unsigned int, 3> > >;
template class lattice_wrapper<float, close_packed_lattice<fixed_vector<float, 2>, fixed_vector<unsigned int, 2> > >;
template class lattice_wrapper<dsfloat, close_packed_lattice<fixed_vector<float, 3>, fixed_vector<unsigned int, 3> > >;
template class lattice_wrapper<dsfloat, close_packed_lattice<fixed_vector<float, 2>, fixed_vector<unsigned int, 2> > >;

} // namespace mdsim
} // namespace gpu
} // namespace positions
} // namespace halmd
