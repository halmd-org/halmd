/*
 * Copyright © 2008-2011  Peter Colberg
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
#include <halmd/mdsim/gpu/velocities/trajectory.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu { namespace velocities
{

template <int dimension, typename float_type>
trajectory<dimension, float_type>::trajectory(
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
void trajectory<dimension, float_type>::set()
{
    for (size_t j = 0, i = 0; j < particle->ntype; i += particle->ntypes[j], ++j) {
        copy(sample->v[j]->begin(), sample->v[j]->end(), &particle->h_v[i]);
    }

#ifdef USE_VERLET_DSFUN
    // erase particle velocity vectors (double-single precision)
    cuda::memset(particle->g_v, 0, particle->g_v.capacity());
#endif

    try {
        cuda::copy(particle->h_v, particle->g_v);
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("[trajectory] failed to copy particle velocities to GPU");
        throw;
    }

    LOG("set particle velocities from trajectory sample");
}

template <int dimension, typename float_type>
void trajectory<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("trajectory_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("gpu")
                [
                    namespace_("velocities")
                    [
                        class_<trajectory, shared_ptr<_Base>, _Base>(class_name.c_str())
                            .def(constructor<
                                 shared_ptr<particle_type>
                               , shared_ptr<sample_type>
                            >())
                    ]
                ]
            ]
        ]
    ];
}


namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
    [
        &trajectory<3, float>::luaopen
    ]
    [
        &trajectory<2, float>::luaopen
    ];
}

} // namespace

// explicit instantiation
template class trajectory<3, float>;
template class trajectory<2, float>;

}}} // namespace mdsim::gpu::velocities

} // namespace halmd
