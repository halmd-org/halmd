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
#include <halmd/observables/host/phase_space.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace host
{

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(
    shared_ptr<sample_type> sample
  , shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
)
  : sample(sample)
  , particle(particle)
  , box(box)
{
}

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
void phase_space<dimension, float_type>::acquire(double time)
{
    if (sample->time == time) {
        LOG_TRACE("[phase_space] sample is up to date");
        return;
    }

    LOG_TRACE("[phase_space] acquire sample");

    for (size_t i = 0; i < particle->nbox; ++i) {
        unsigned int type = particle->type[i];
        unsigned int tag = particle->tag[i];

        // periodically extended particle position
        assert(type < sample->r.size());
        assert(tag < sample->r[type]->size());
        (*sample->r[type])[tag] = particle->r[i] + element_prod(particle->image[i], box->length());

        // particle velocity
        assert(type < sample->v.size());
        assert(tag < sample->v[type]->size());
        (*sample->v[type])[tag] = particle->v[i];
    }
    sample->time = time;
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("phase_space_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<phase_space, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .def(constructor<
                         shared_ptr<sample_type>
                       , shared_ptr<particle_type>
                       , shared_ptr<box_type>
                    >())
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_phase_space(lua_State* L)
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

}} // namespace observables::host

} // namespace halmd
