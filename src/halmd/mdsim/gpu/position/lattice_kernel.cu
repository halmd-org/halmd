/*
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

#include <boost/mpl/int.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/mdsim/gpu/position/lattice_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>
#include <halmd/utility/gpu/thread.cuh>

using namespace boost;
using namespace halmd::numeric::gpu;
using namespace halmd::mdsim::gpu::particle_kernel;

namespace halmd
{
namespace mdsim { namespace gpu { namespace position
{
namespace lattice_kernel
{

using boost::mpl::int_;

/**
 * place particles on a face centered cubic lattice (fcc)
 */
template <typename vector_type>
__device__ typename enable_if<is_same<int_<3>, int_<vector_type::static_size> > >::type
fcc(vector_type& r, uint n)
{
    // compose primitive vectors from 1-dimensional index
    r.x = ((GTID >> 2) % n) + ((GTID ^ (GTID >> 1)) & 1) / 2.f;
    r.y = ((GTID >> 2) / n % n) + (GTID & 1) / 2.f;
    r.z = ((GTID >> 2) / n / n) + (GTID & 2) / 4.f;
}

template <typename vector_type>
__device__ typename enable_if<is_same<int_<2>, int_<vector_type::static_size> > >::type
fcc(vector_type& r, uint n)
{
    r.x = ((GTID >> 1) % n) + (GTID & 1) / 2.f;
    r.y = ((GTID >> 1) / n) + (GTID & 1) / 2.f;
}

/**
 * place particles on a simple cubic lattice (sc)
 */
template <typename vector_type>
__device__ typename enable_if<is_same<int_<3>, int_<vector_type::static_size> > >::type
sc(vector_type& r, uint n)
{
    r.x = (GTID % n) + 0.5f;
    r.y = (GTID / n % n) + 0.5f;
    r.z = (GTID / n / n) + 0.5f;
}

template <typename vector_type>
__device__ typename enable_if<is_same<int_<2>, int_<vector_type::static_size> > >::type
sc(vector_type& r, uint n)
{
    r.x = (GTID % n) + 0.5f;
    r.y = (GTID / n) + 0.5f;
}

template <
    typename vector_type
  , void (*primitive)(vector_type&, uint)
>
__global__ void lattice(float4* g_r, uint n, float box)
{
    vector_type r;
    unsigned int type;
    tie(r, type) = untagged<vector_type>(g_r[GTID]);

    // compute primitive lattice vector
    primitive(r, n);
    r *= box / n;

    g_r[GTID] = tagged(r, type);
}

} // namespace lattice_kernel

template <int dimension>
lattice_wrapper<dimension> const lattice_wrapper<dimension>::kernel = {
    lattice_kernel::lattice<blas::vector<float, dimension>, lattice_kernel::fcc>
  , lattice_kernel::lattice<blas::vector<float, dimension>, lattice_kernel::sc>
};

template class lattice_wrapper<3>;
template class lattice_wrapper<2>;

}}} // namespace mdsim::gpu::position

} // namespace halmd
