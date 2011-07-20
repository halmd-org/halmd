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

#include <halmd/mdsim/host/positions/phase_space.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace positions {

using namespace boost;
using namespace std;

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type const> box
  , shared_ptr<sample_type const> sample
  , shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , sample_(sample)
  , logger_(logger)
{
}

/**
 * set particle positions
 */
template <int dimension, typename float_type>
void phase_space<dimension, float_type>::set()
{
    LOG("set particle positions from phase space sample");

    // assign particle coordinates
    for (size_t j = 0, i = 0; j < particle_->ntype; i += particle_->ntypes[j], ++j) {
        assert(sample_->r[j]->size() + i <= particle_->r.size());
        copy(sample_->r[j]->begin(), sample_->r[j]->end(), particle_->r.begin() + i);
    }

    // shift particle positions to range (-L/2, L/2)
    for (size_t i = 0; i < particle_->nbox; ++i) {
        vector_type& r = particle_->r[i];
        vector_type& image = particle_->image[i] = 0;

        // The host implementation of reduce_periodic wraps the position at
        // most once around the box. This is more efficient during the
        // simulation. For setting arbitrary particle positions, however, we
        // must ensure that the final position is actually inside the periodic
        // box.
        vector_type a;
        do {
            image += (a = box_->reduce_periodic(r));
        } while(a != vector_type(0));
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
                namespace_("positions")
                [
                    class_<phase_space, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                             shared_ptr<particle_type>
                           , shared_ptr<box_type const>
                           , shared_ptr<sample_type const>
                           , shared_ptr<logger_type>
                        >())
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_positions_phase_space(lua_State* L)
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

} // namespace mdsim
} // namespace host
} // namespace positions
} // namespace halmd
