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

#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {

template <typename thermodynamics_type>
static function<double ()>
wrap_nparticle(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::nparticle, thermodynamics);
}

template <typename thermodynamics_type>
static function<double ()>
wrap_density(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::density, thermodynamics);
}

template <typename thermodynamics_type>
static function<double ()>
wrap_en_tot(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::en_tot, thermodynamics);
}

template <typename thermodynamics_type>
static function<double ()>
wrap_en_pot(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::en_pot, thermodynamics);
}

template <typename thermodynamics_type>
static function<double ()>
wrap_en_kin(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::en_kin, thermodynamics);
}

template <typename thermodynamics_type>
static function<typename thermodynamics_type::vector_type ()>
wrap_v_cm(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::v_cm, thermodynamics);
}

template <typename thermodynamics_type>
static function<double ()>
wrap_temp(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::temp, thermodynamics);
}

template <typename thermodynamics_type>
static function<double ()>
wrap_pressure(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::pressure, thermodynamics);
}

template <typename thermodynamics_type>
static function<double ()>
wrap_virial(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::virial, thermodynamics);
}

template <typename thermodynamics_type>
static function<double ()>
wrap_hypervirial(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::hypervirial, thermodynamics);
}

template <int dimension>
void thermodynamics<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("thermodynamics_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<thermodynamics, shared_ptr<thermodynamics> >(class_name.c_str())
                .property("nparticle", &wrap_nparticle<thermodynamics>)
                .property("density", &wrap_density<thermodynamics>)
                .property("en_kin", &wrap_en_kin<thermodynamics>)
                .property("en_pot", &wrap_en_pot<thermodynamics>)
                .property("en_tot", &wrap_en_tot<thermodynamics>)
                .property("pressure", &wrap_pressure<thermodynamics>)
                .property("temp", &wrap_temp<thermodynamics>)
                .property("v_cm", &wrap_v_cm<thermodynamics>)
                .property("virial", &wrap_virial<thermodynamics>)
                .property("hypervirial", &wrap_hypervirial<thermodynamics>)
                .def("clear_cache", &thermodynamics::clear_cache)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_thermodynamics(lua_State* L)
{
    thermodynamics<3>::luaopen(L);
    thermodynamics<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class thermodynamics<3>;
template class thermodynamics<2>;

} // namespace observables
} // namespace halmd
