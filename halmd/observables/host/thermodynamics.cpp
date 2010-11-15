/*
 * Copyright © 2010  Felix Höfling
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
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace host
{

template <int dimension, typename float_type>
thermodynamics<dimension, float_type>::thermodynamics(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<force_type> force
)
  : _Base(box)
  // dependency injection
  , particle(particle)
  , force(force)
{
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_kin() const
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
typename thermodynamics<dimension, float_type>::vector_type thermodynamics<dimension, float_type>::v_cm() const
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
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static string class_name("thermodynamics_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("observables")
            [
                namespace_("host")
                [
                    class_<thermodynamics, shared_ptr<_Base_Base>, bases<_Base, _Base_Base> >(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type>
                          , shared_ptr<box_type>
                          , shared_ptr<force_type>
                        >())
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &thermodynamics<3, double>::luaopen
    ]
    [
        &thermodynamics<2, double>::luaopen
    ];
#else
    [
        &thermodynamics<3, float>::luaopen
    ]
    [
        &thermodynamics<2, float>::luaopen
    ];
#endif
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class thermodynamics<3, double>;
template class thermodynamics<2, double>;
#else
template class thermodynamics<3, float>;
template class thermodynamics<2, float>;
#endif

}} // namespace observables::host

} // namespace halmd
