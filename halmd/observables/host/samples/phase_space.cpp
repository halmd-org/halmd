/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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
#include <boost/lexical_cast.hpp>
#include <limits>
#include <stdexcept>
#include <string>

#include <halmd/observables/host/samples/phase_space.hpp>
#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {
namespace samples {

template <int dimension, typename float_type>
static int wrap_dimension(phase_space<dimension, float_type> const&)
{
    return dimension;
}

template <int dimension, typename float_type>
typename phase_space<dimension, float_type>::sample_vector const&
phase_space<dimension, float_type>::position() const
{
    return *r;
}

template <int dimension, typename float_type>
typename phase_space<dimension, float_type>::sample_vector const&
phase_space<dimension, float_type>::velocity() const
{
    return *v;
}

template <int dimension, typename float_type>
typename phase_space<dimension, float_type>::type_vector const&
phase_space<dimension, float_type>::types() const
{
    return *type;
}

template <int dimension, typename float_type>
typename phase_space<dimension, float_type>::sample_vector&
phase_space<dimension, float_type>::position()
{
    return *r;
}

template <int dimension, typename float_type>
typename phase_space<dimension, float_type>::sample_vector&
phase_space<dimension, float_type>::velocity()
{
    return *v;
}

template <int dimension, typename float_type>
typename phase_space<dimension, float_type>::type_vector&
phase_space<dimension, float_type>::types()
{
    return *type;
}

template <typename phase_space_type>
static function<typename phase_space_type::sample_vector& ()>
wrap_position(shared_ptr<phase_space_type> phase_space)
{
    typedef typename phase_space_type::sample_vector& (phase_space_type::*getter_type)();
    return bind(static_cast<getter_type>(&phase_space_type::position), phase_space);
}

template <typename phase_space_type>
static function<typename phase_space_type::sample_vector& ()>
wrap_velocity(shared_ptr<phase_space_type> phase_space)
{
    typedef typename phase_space_type::sample_vector& (phase_space_type::*getter_type)();
    return bind(static_cast<getter_type>(&phase_space_type::velocity), phase_space);
}

template <typename phase_space_type>
static function<typename phase_space_type::type_vector& ()>
wrap_types(shared_ptr<phase_space_type> phase_space)
{
    typedef typename phase_space_type::type_vector& (phase_space_type::*getter_type)();
    return bind(static_cast<getter_type>(&phase_space_type::types), phase_space);
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name("phase_space_" + lexical_cast<string>(dimension) + "_" + demangled_name<float_type>() + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                namespace_("samples")
                [
                    class_<phase_space, shared_ptr<phase_space> >(class_name.c_str())
                        .def(constructor<unsigned int>())
                        .property("dimension", &wrap_dimension<dimension, float_type>)
                        .property("position", &wrap_position<phase_space>)
                        .property("velocity", &wrap_velocity<phase_space>)
                        .property("types", &wrap_types<phase_space>)
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
    observables::samples::blocking_scheme<phase_space<3, double> >::luaopen(L);
    observables::samples::blocking_scheme<phase_space<2, double> >::luaopen(L);
#endif
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
    observables::samples::blocking_scheme<phase_space<3, float> >::luaopen(L);
    observables::samples::blocking_scheme<phase_space<2, float> >::luaopen(L);
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class phase_space<3, double>;
template class phase_space<2, double>;
#endif
template class phase_space<3, float>;
template class phase_space<2, float>;

} // namespace samples
} // namespace host

namespace samples
{

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class blocking_scheme<host::samples::phase_space<3, double> >;
template class blocking_scheme<host::samples::phase_space<2, double> >;
#endif
template class blocking_scheme<host::samples::phase_space<3, float> >;
template class blocking_scheme<host::samples::phase_space<2, float> >;

} // namespace samples
} // namespace observables
} // namespace halmd
