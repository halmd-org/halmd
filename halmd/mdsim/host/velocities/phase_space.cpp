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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/velocities/phase_space.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace velocities
{

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(
    shared_ptr<particle_type> particle
  , shared_ptr<sample_type> sample
)
  : _Base(particle)
  // dependency injection
  , particle(particle)
  , sample(sample)
{
}

/**
 * set particle velocities
 */
template <int dimension, typename float_type>
void phase_space<dimension, float_type>::set()
{
    LOG("set particle velocities from phase space sample");

    for (size_t j = 0, i = 0; j < particle->ntype; i += particle->ntypes[j], ++j) {
        assert(sample->v[j]->size() + i <= particle->v.size());
        copy(sample->v[j]->begin(), sample->v[j]->end(), particle->v.begin() + i);
    }
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("phase_space_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("velocities")
                [
                    class_<phase_space, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                             shared_ptr<particle_type>
                           , shared_ptr<sample_type>
                        >())
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_velocities_phase_space(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    phase_space<3, double>::luaopen(L);
    phase_space<2, double>::luaopen(L);
#else
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class phase_space<3, double>;
template class phase_space<2, double>;
#else
template class phase_space<3, float>;
template class phase_space<2, float>;
#endif

}}} // namespace mdsim::host::velocities

} // namespace halmd