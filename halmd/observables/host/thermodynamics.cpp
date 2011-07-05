/*
 * Copyright © 2010-2011  Felix Höfling
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

#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
thermodynamics<dimension, float_type>::thermodynamics(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<clock_type> clock
  , shared_ptr<force_type> force
)
  : _Base(box, clock)
  // dependency injection
  , particle(particle)
  , force(force)
{
}

/**
 * preparations before computation of forces
 *
 * set flag of force module to compute auxiliary
 * variables like potential energy, stress tensor,
 * and hypervirial
 */
template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::prepare()
{
    force->aux_enable();
}

/**
 * call sample() from base class and
 * unset flag for auxiliary variables of force module at the end
 */
template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::sample(uint64_t step)
{
    _Base::sample(step);
    force->aux_disable();
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_kin()
{
    // compute mean-square velocity
    double vv = 0;
    BOOST_FOREACH(vector_type const& v, particle->v) {
        // assuming unit mass for all particle types
        vv += inner_prod(v, v);
    }
    return .5 * vv / particle->nbox;
}

template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::vector_type thermodynamics<dimension, float_type>::v_cm()
{
    // compute mean velocity
    vector_type v_cm_(0.);
    BOOST_FOREACH(vector_type const& v, particle->v) {
        v_cm_ += v;
    }
    return v_cm_ / particle->nbox;
}

template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("thermodynamics_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<thermodynamics, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .def(constructor<
                        shared_ptr<particle_type>
                      , shared_ptr<box_type>
                      , shared_ptr<clock_type>
                      , shared_ptr<force_type>
                    >())
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_thermodynamics(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    thermodynamics<3, double>::luaopen(L);
    thermodynamics<2, double>::luaopen(L);
#else
    thermodynamics<3, float>::luaopen(L);
    thermodynamics<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class thermodynamics<3, double>;
template class thermodynamics<2, double>;
#else
template class thermodynamics<3, float>;
template class thermodynamics<2, float>;
#endif

} // namespace observables
} // namespace host
} // namespace halmd
