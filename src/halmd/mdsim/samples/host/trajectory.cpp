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

#include <halmd/mdsim/samples/host/trajectory.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace samples { namespace host
{

template <int dimension, typename float_type>
trajectory<dimension, float_type>::trajectory(
    shared_ptr<particle_type> particle
)
  // dependency injection
  : particle(particle)
  // allocate sample pointers
  , r(particle->ntype)
  , v(particle->ntype)
{
    for (size_t i = 0; i < particle->ntype; ++i) {
        r[i].reset(new sample_vector(particle->ntypes[i]));
        v[i].reset(new sample_vector(particle->ntypes[i]));
    }
}

template <typename T>
static void register_lua(char const* class_name)
{
    typedef typename T::particle_type particle_type;

    using namespace luabind;
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("samples")
                [
                    namespace_("host")
                    [
                        class_<T, shared_ptr<T> >(class_name)
                            .def("acquire", &T::acquire)
                    ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
#ifndef USE_HOST_SINGLE_PRECISION
    register_lua<trajectory<3, double> >("trajectory_3_double_");
    register_lua<trajectory<2, double> >("trajectory_2_double_");
#endif
    register_lua<trajectory<3, float> >("trajectory_3_float_");
    register_lua<trajectory<2, float> >("trajectory_2_float_");
}

#ifndef USE_HOST_SINGLE_PRECISION
template class trajectory<3, double>;
template class trajectory<2, double>;
#endif
template class trajectory<3, float>;
template class trajectory<2, float>;

}}} // namespace mdsim::samples::host

} // namespace halmd
