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
#include <cmath>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/integrators/verlet.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace integrators
{

template <int dimension, typename float_type>
verlet<dimension, float_type>::verlet(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , double timestep
)
  // dependency injection
  : particle(particle)
  , box(box)
{
    this->timestep(timestep);
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::timestep(double timestep)
{
  timestep_ = timestep;
  timestep_half_ = 0.5 * timestep;

  LOG("integration timestep: " << timestep_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::integrate()
{
    for (size_t i = 0; i < particle->nbox; ++i) {
        vector_type& v = particle->v[i] += particle->f[i] * timestep_half_;
        vector_type& r = particle->r[i] += v * timestep_;
        // enforce periodic boundary conditions
        // TODO: reduction is now to (-L/2, L/2) instead of (0, L) as before
        // check that this is OK
        particle->image[i] += box->reduce_periodic(r);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::finalize()
{
    for (size_t i = 0; i < particle->nbox; ++i) {
        particle->v[i] += particle->f[i] * timestep_half_;
    }
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(verlet<dimension, float_type> const&)
{
    return verlet<dimension, float_type>::module_name();
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("integrators")
                [
                    class_<verlet, shared_ptr<_Base>, bases<_Base> >(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type>
                          , shared_ptr<box_type>
                          , double>()
                        )
                        .property("module_name", &module_name_wrapper<dimension, float_type>)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_integrators_verlet(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    verlet<3, double>::luaopen(L);
    verlet<2, double>::luaopen(L);
#else
    verlet<3, float>::luaopen(L);
    verlet<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class verlet<3, double>;
template class verlet<2, double>;
#else
template class verlet<3, float>;
template class verlet<2, float>;
#endif

}}} // namespace mdsim::host::integrators

} // namespace halmd
