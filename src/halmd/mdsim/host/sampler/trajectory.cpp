/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/sampler/trajectory.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace sampler
{

template <int dimension, typename float_type>
trajectory<dimension, float_type>::trajectory(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<core_type> core
)
  : _Base(particle)     //< mdsim::particle
  , particle(particle)  //< mdsim::host:particle
  , box(box)
  , core(core)
{
}

/**
 * Sample trajectory
 */
template <int dimension, typename float_type>
void trajectory<dimension, float_type>::acquire()
{
    for (size_t i = 0; i < particle->nbox; ++i) {
        // periodically extended particle position
        (*r[particle->type[i]])[particle->tag[i]] = particle->r[i] + element_prod(particle->image[i], box->length());
        // particle velocity
        (*v[particle->type[i]])[particle->tag[i]] = particle->v[i];
    }
    time = core->time();
}

template <typename T>
static void register_lua(char const* class_name)
{
    typedef typename T::_Base _Base;
    typedef typename T::particle_type particle_type;
    typedef typename T::box_type box_type;
    typedef typename T::core_type core_type;

    using namespace luabind;
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    namespace_("sampler")
                    [
                        class_<T, shared_ptr<_Base>, _Base>(class_name)
                            .def(constructor<
                                 shared_ptr<particle_type>
                               , shared_ptr<box_type>
                               , shared_ptr<core_type>
                            >())
                    ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
#ifndef USE_HOST_SINGLE_PRECISION
    register_lua<trajectory<3, double> >("trajectory_3_");
    register_lua<trajectory<2, double> >("trajectory_2_");
#else
    register_lua<trajectory<3, float> >("trajectory_3_");
    register_lua<trajectory<2, float> >("trajectory_2_");
#endif
}

}}} // namespace mdsim::host::sampler

} // namespace halmd
