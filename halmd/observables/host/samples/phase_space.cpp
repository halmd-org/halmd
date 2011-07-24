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
#include <boost/bind.hpp>
#include <limits>
#include <stdexcept>

#include <halmd/observables/host/samples/phase_space.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {
namespace samples {

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(vector<unsigned int> ntypes)
  // allocate sample pointers
  : r(ntypes.size())
  , v(ntypes.size())
  // initialise attributes
  , step(numeric_limits<uint64_t>::max())
{
    for (size_t i = 0; i < ntypes.size(); ++i) {
        r[i].reset(new sample_vector(ntypes[i]));
        v[i].reset(new sample_vector(ntypes[i]));
    }
}

template <int dimension, typename float_type>
static int wrap_dimension(phase_space<dimension, float_type> const&)
{
    return dimension;
}

template <int dimension, typename float_type>
typename phase_space<dimension, float_type>::sample_vector const&
phase_space<dimension, float_type>::position(unsigned int type) const
{
    if (!(type < r.size())) {
        throw invalid_argument("particle type");
    }
    return *r[type];
}

template <int dimension, typename float_type>
typename phase_space<dimension, float_type>::sample_vector const&
phase_space<dimension, float_type>::velocity(unsigned int type) const
{
    if (!(type < v.size())) {
        throw invalid_argument("particle type");
    }
    return *v[type];
}

template <int dimension, typename float_type>
typename phase_space<dimension, float_type>::sample_vector&
phase_space<dimension, float_type>::position(unsigned int type)
{
    if (!(type < r.size())) {
        throw invalid_argument("particle type");
    }
    return *r[type];
}

template <int dimension, typename float_type>
typename phase_space<dimension, float_type>::sample_vector&
phase_space<dimension, float_type>::velocity(unsigned int type)
{
    if (!(type < v.size())) {
        throw invalid_argument("particle type");
    }
    return *v[type];
}

template <typename phase_space_type>
static function<typename phase_space_type::sample_vector& ()>
wrap_position(shared_ptr<phase_space_type> phase_space, unsigned int type)
{
    typedef typename phase_space_type::sample_vector& (phase_space_type::*getter_type)(unsigned int);
    return bind(static_cast<getter_type>(&phase_space_type::position), phase_space, type);
}

template <typename phase_space_type>
static function<typename phase_space_type::sample_vector& ()>
wrap_velocity(shared_ptr<phase_space_type> phase_space, unsigned int type)
{
    typedef typename phase_space_type::sample_vector& (phase_space_type::*getter_type)(unsigned int);
    return bind(static_cast<getter_type>(&phase_space_type::velocity), phase_space, type);
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("phase_space_" + lexical_cast<string>(dimension) + "_" + demangled_name<float_type>() + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                namespace_("samples")
                [
                    class_<phase_space, shared_ptr<phase_space> >(class_name.c_str())
                        .def(constructor<vector<unsigned int> >())
                        .property("dimension", &wrap_dimension<dimension, float_type>)
                        .def("position", &wrap_position<phase_space>)
                        .def("velocity", &wrap_velocity<phase_space>)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_samples_phase_space(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    phase_space<3, double>::luaopen(L);
    phase_space<2, double>::luaopen(L);
#endif
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
    return 0;
}

#ifndef USE_HOST_SINGLE_PRECISION
template class phase_space<3, double>;
template class phase_space<2, double>;
#endif
template class phase_space<3, float>;
template class phase_space<2, float>;

} // namespace observables
} // namespace host
} // namespace samples
} // namespace halmd
