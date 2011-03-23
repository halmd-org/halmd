/*
 * Copyright © 2008-2009  Peter Colberg and Felix Höfling
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

#include <boost/mpl/int.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/mdsim/gpu/positions/lattice_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/gpu/variant.cuh>

using namespace boost;
using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::utility::gpu;

namespace halmd
{
namespace mdsim { namespace gpu { namespace positions
{
namespace lattice_kernel
{

using boost::mpl::int_;

/** edge lengths of cuboid slab */
static __constant__ variant<map<pair<int_<3>, float3>, pair<int_<2>, float2> > > offset_;
/** number of cells per dimension */
static __constant__ variant<map<pair<int_<3>, uint3>, pair<int_<2>, uint2> > > ncell_;

/**
 * place particles on a face centered cubic lattice (fcc)
 */
template <typename vector_type, typename index_type>
__device__ typename enable_if<is_same<int_<3>, int_<vector_type::static_size> > >::type
fcc(unsigned int i, index_type const& n, vector_type& r)
{
    // compose primitive vectors from 1-dimensional index
    r[0] = ((i >> 2) % n[0]) + ((i ^ (i >> 1)) & 1) / 2.f;
    r[1] = ((i >> 2) / n[0] % n[1]) + (i & 1) / 2.f;
    r[2] = ((i >> 2) / n[0] / n[1]) + (i & 2) / 4.f;
}

template <typename vector_type, typename index_type>
__device__ typename enable_if<is_same<int_<2>, int_<vector_type::static_size> > >::type
fcc(unsigned int i, index_type const& n, vector_type& r)
{
    r[0] = ((i >> 1) % n[0]) + (i & 1) / 2.f;
    r[1] = ((i >> 1) / n[0]) + (i & 1) / 2.f;
}

/**
 * place particles on a simple cubic lattice (sc)
 */
template <typename vector_type, typename index_type>
__device__ typename enable_if<is_same<int_<3>, int_<vector_type::static_size> > >::type
sc(unsigned int i, index_type const& n, vector_type& r)
{
    r[0] = (i % n[0]) + 0.5f;
    r[1] = (i / n[0] % n[1]) + 0.5f;
    r[2] = (i / n[0] / n[1]) + 0.5f;
}

template <typename vector_type, typename index_type>
__device__ typename enable_if<is_same<int_<2>, int_<vector_type::static_size> > >::type
sc(unsigned int i, index_type const& n, vector_type& r)
{
    r[0] = (i % n[0]) + 0.5f;
    r[1] = (i / n[0]) + 0.5f;
}

template <
    typename vector_type
  , void (*primitive)(
        unsigned int
      , fixed_vector<unsigned int, vector_type::static_size> const&
      , vector_type&
    )
>
__global__ void lattice(float4* g_r, uint npart, float a, uint skip)
{
    enum { dimension = vector_type::static_size };

    for (uint i = GTID; i < npart; i += GTDIM) {

        // load particle type
        vector_type r;
        unsigned int type;
        tie(r, type) = untagged<vector_type>(g_r[i]);

        // introduce a vacancy after every (skip - 1) particles
        uint nvacancies = (skip > 1) ? (i / (skip - 1)) : 0;

        // compute primitive lattice vector
        primitive(i + nvacancies, get<dimension>(ncell_), r);

        // scale with lattice constant
        r *= a;

        // shift origin of lattice to offset
        r += vector_type(get<dimension>(offset_));

        g_r[i] = tagged(r, type);
    }
}

} // namespace lattice_kernel

template <int dimension>
lattice_wrapper<dimension> const lattice_wrapper<dimension>::kernel = {
    get<dimension>(lattice_kernel::offset_)
  , get<dimension>(lattice_kernel::ncell_)
  , lattice_kernel::lattice<fixed_vector<float, dimension>, lattice_kernel::fcc>
  , lattice_kernel::lattice<fixed_vector<float, dimension>, lattice_kernel::sc>
};

template class lattice_wrapper<3>;
template class lattice_wrapper<2>;

}}} // namespace mdsim::gpu::positions

} // namespace halmd
