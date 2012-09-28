/*
 * Copyright © 2012 Peter Colberg
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

namespace halmd {
namespace observables {

template <typename thermodynamics_type>
static std::function<double ()>
wrap_nparticle(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->nparticle();
    };
}

template <typename thermodynamics_type>
static std::function<double ()>
wrap_density(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->density();
    };
}

template <typename thermodynamics_type>
static std::function<double ()>
wrap_en_tot(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->en_tot();
    };
}

template <typename thermodynamics_type>
static std::function<double ()>
wrap_en_pot(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->en_pot();
    };
}

template <typename thermodynamics_type>
static std::function<double ()>
wrap_en_kin(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->en_kin();
    };
}

template <typename thermodynamics_type>
static std::function<typename thermodynamics_type::vector_type ()>
wrap_v_cm(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->v_cm();
    };
}

template <typename thermodynamics_type>
static std::function<typename thermodynamics_type::vector_type ()>
wrap_r_cm(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->r_cm();
    };
}

template <typename thermodynamics_type>
static std::function<double ()>
wrap_mean_mass(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->mean_mass();
    };
}

template <typename thermodynamics_type>
static std::function<double ()>
wrap_temp(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->temp();
    };
}

template <typename thermodynamics_type>
static std::function<double ()>
wrap_pressure(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->pressure();
    };
}

template <typename thermodynamics_type>
static std::function<double ()>
wrap_virial(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->virial();
    };
}

template <typename thermodynamics_type>
static std::function<double ()>
wrap_hypervirial(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->hypervirial();
    };
}

template <int dimension>
void thermodynamics<dimension>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L)
    [
        class_<thermodynamics>()
            .property("nparticle", &wrap_nparticle<thermodynamics>)
            .property("density", &wrap_density<thermodynamics>)
            .property("en_kin", &wrap_en_kin<thermodynamics>)
            .property("en_pot", &wrap_en_pot<thermodynamics>)
            .property("en_tot", &wrap_en_tot<thermodynamics>)
            .property("pressure", &wrap_pressure<thermodynamics>)
            .property("temp", &wrap_temp<thermodynamics>)
            .property("v_cm", &wrap_v_cm<thermodynamics>)
            .property("r_cm", &wrap_r_cm<thermodynamics>)
            .property("mean_mass", &wrap_mean_mass<thermodynamics>)
            .property("virial", &wrap_virial<thermodynamics>)
            .property("hypervirial", &wrap_hypervirial<thermodynamics>)
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
