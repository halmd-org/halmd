/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
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

template <typename sample_type>
static typename sample_type::position_array_type const&
position(function<shared_ptr<sample_type const> ()> const& slot, shared_ptr<sample_type const>& sample)
{
    sample = slot();
    return sample->position();
}

template <typename sample_type>
static function<typename sample_type::position_array_type const& ()>
wrap_position(function<shared_ptr<sample_type const> ()> const& slot)
{
    return bind(&position<sample_type>, slot, shared_ptr<sample_type const>());
}

template <typename sample_type>
static typename sample_type::velocity_array_type const&
velocity(function<shared_ptr<sample_type const> ()> const& slot, shared_ptr<sample_type const>& sample)
{
    sample = slot();
    return sample->velocity();
}

template <typename sample_type>
static function<typename sample_type::velocity_array_type const& ()>
wrap_velocity(function<shared_ptr<sample_type const> ()> const& slot)
{
    return bind(&velocity<sample_type>, slot, shared_ptr<sample_type const>());
}

template <typename sample_type>
static typename sample_type::species_array_type const&
species(function<shared_ptr<sample_type const> ()> const& slot, shared_ptr<sample_type const>& sample)
{
    sample = slot();
    return sample->species();
}

template <typename sample_type>
static function<typename sample_type::species_array_type const& ()>
wrap_species(function<shared_ptr<sample_type const> ()> const& slot)
{
    return bind(&species<sample_type>, slot, shared_ptr<sample_type const>());
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name("phase_space_" + lexical_cast<string>(dimension) + "_" + demangled_name<float_type>());
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                namespace_("samples")
                [
                    class_<phase_space>(class_name.c_str())
                        .property("step", &phase_space::step)
                ]
            ]
          , namespace_("samples")
            [
                namespace_("phase_space")
                [
                    def("position", &wrap_position<phase_space>)
                  , def("velocity", &wrap_velocity<phase_space>)
                  , def("species", &wrap_species<phase_space>)
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
