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

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <exception>
#include <numeric>

#include <halmd/algorithm/host/permute.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host
{

/**
 * Allocate microscopic system state.
 *
 * @param particles number of particles per type or species
 */
template <unsigned int dimension, typename float_type>
particle<dimension, float_type>::particle(vector<unsigned int> const& particles)
  : _Base(particles)
  // allocate particle storage
  , r(nbox)
  , image(nbox)
  , v(nbox)
  , f(nbox)
  , tag(nbox)
  , type(nbox)
  , neighbour(nbox)
{
}

/**
 * set particle tags and types
 */
template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::set()
{
    // handle each type separately
    for (size_t i = 0, j = 0; j < ntype; i += ntypes[j], ++j) {
        // assign particle types
        fill_n(type.begin() + i, ntypes[j], j);
        // assign particle tags
        copy(
            counting_iterator<size_t>(0)
          , counting_iterator<size_t>(ntypes[j])
          , tag.begin() + i
        );
    }
}

/**
 * Rearrange particles in memory according to an integer index sequence
 *
 * The neighbour lists must be rebuilt after calling this function!
 */
template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::rearrange(std::vector<unsigned int> const& index)
{
    algorithm::host::permute(r.begin(), r.end(), index.begin());
    algorithm::host::permute(image.begin(), image.end(), index.begin());
    algorithm::host::permute(v.begin(), v.end(), index.begin());
    // no permutation of forces
    algorithm::host::permute(tag.begin(), tag.end(), index.begin());
    algorithm::host::permute(type.begin(), type.end(), index.begin());
    // no permutation of neighbour lists
}

template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("particle_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    class_<particle, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<vector<unsigned int> const&>())
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_particle(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    particle<3, double>::luaopen(L);
    particle<2, double>::luaopen(L);
#else
    particle<3, float>::luaopen(L);
    particle<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class particle<3, double>;
template class particle<2, double>;
#else
template class particle<3, float>;
template class particle<2, float>;
#endif

}} // namespace mdsim::host

} // namespace halmd
