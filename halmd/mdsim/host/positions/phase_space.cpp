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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/positions/phase_space.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace positions
{

using namespace boost;
using namespace std;

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<sample_type> sample
)
  // dependency injection
  : particle(particle)
  , box(box)
  , sample(sample)
{
}

/**
 * set particle positions
 */
template <int dimension, typename float_type>
void phase_space<dimension, float_type>::set()
{
    // assign particle coordinates
    for (size_t j = 0, i = 0; j < particle->ntype; i += particle->ntypes[j], ++j) {
        copy(sample->r[j]->begin(), sample->r[j]->end(), &particle->r[i]);
    }

    // shift particle positions to range (-L/2, L/2)
    for (size_t i = 0; i < particle->nbox; ++i) {
        box->reduce_periodic(particle->r[i]);
    }

    // assign particle image vectors
    fill(particle->image.begin(), particle->image.end(), 0);

    LOG("set particle positions from phase space sample");
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("phase_space_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    namespace_("positions")
                    [
                        class_<phase_space, shared_ptr<_Base>, _Base>(class_name.c_str())
                            .def(constructor<
                                 shared_ptr<particle_type>
                               , shared_ptr<box_type>
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
    lua_wrapper::register_(1) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &phase_space<3, double>::luaopen
    ]
    [
        &phase_space<2, double>::luaopen
    ];
#else
    [
        &phase_space<3, float>::luaopen
    ]
    [
        &phase_space<2, float>::luaopen
    ];
#endif
}

} // namespace

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class phase_space<3, double>;
template class phase_space<2, double>;
#else
template class phase_space<3, float>;
template class phase_space<2, float>;
#endif

}}} // namespace mdsim::host::positions

} // namespace halmd
