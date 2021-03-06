/*
 * Copyright © 2016 Felix Höfling
 * Copyright © 2012 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {

template <typename thermodynamics_type>
static std::function<unsigned int ()>
wrap_particle_number(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->particle_number();
    };
}

template <typename thermodynamics_type>
static std::function<double ()>
wrap_volume(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->volume();
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
wrap_total_force(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->total_force();
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
static std::function<typename thermodynamics_type::stress_tensor_type ()>
wrap_stress_tensor(std::shared_ptr<thermodynamics_type> self)
{
    return [=]() {
        return self->stress_tensor();
    };
}

template <int dimension>
void thermodynamics<dimension>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L)
    [
        class_<thermodynamics>()
            .property("particle_number", &wrap_particle_number<thermodynamics>)
            .property("volume", &wrap_volume<thermodynamics>)
            .property("density", &wrap_density<thermodynamics>)
            .property("kinetic_energy", &wrap_en_kin<thermodynamics>)
            .property("potential_energy", &wrap_en_pot<thermodynamics>)
            .property("internal_energy", &wrap_en_tot<thermodynamics>)
            .property("pressure", &wrap_pressure<thermodynamics>)
            .property("temperature", &wrap_temp<thermodynamics>)
            .property("total_force", &wrap_total_force<thermodynamics>)
            .property("center_of_mass_velocity", &wrap_v_cm<thermodynamics>)
            .property("center_of_mass", &wrap_r_cm<thermodynamics>)
            .property("mean_mass", &wrap_mean_mass<thermodynamics>)
            .property("virial", &wrap_virial<thermodynamics>)
            .property("stress_tensor", &wrap_stress_tensor<thermodynamics>)
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
