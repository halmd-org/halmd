/*
 * Copyright © 2010-2016 Felix Höfling
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2010-2012 Peter Colberg
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

#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
thermodynamics<dimension, float_type>::thermodynamics(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<particle_group_type> group
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , group_(group)
  , box_(box)
  , logger_(logger)
{
}

template <int dimension, typename float_type>
unsigned int thermodynamics<dimension, float_type>::particle_number() const
{
    return *group_->size();
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::volume() const
{
    return box_->volume();
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_kin()
{
    cache<velocity_array_type> const& velocity_cache = particle_->velocity();
    cache<mass_array_type> const& mass_cache = particle_->mass();
    cache<size_type> const& group_cache = group_->size();

    if (en_kin_cache_ != std::tie(velocity_cache, mass_cache, group_cache)) {
        LOG_TRACE("acquire kinetic energy");
        scoped_timer_type timer(runtime_.en_kin);
        en_kin_ = get_mean_en_kin(*particle_, *group_);
        en_kin_cache_ = std::tie(velocity_cache, mass_cache, group_cache);
    }
    return en_kin_;
}

template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::vector_type const&
thermodynamics<dimension, float_type>::total_force()
{
    cache<velocity_array_type> const& force_cache = particle_->force();
    cache<size_type> const& group_cache = group_->size();

    if (force_cache_ != std::tie(force_cache, group_cache)) {
        LOG_TRACE("acquire total force");
        scoped_timer_type timer(runtime_.force);
        force_ = get_total_force(*particle_, *group_);
        force_cache_ = std::tie(force_cache, group_cache);
    }
    return force_;
}

template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::vector_type const&
thermodynamics<dimension, float_type>::v_cm()
{
    cache<velocity_array_type> const& velocity_cache = particle_->velocity();
    cache<mass_array_type> const& mass_cache = particle_->mass();
    cache<size_type> const& group_cache = group_->size();

    if (v_cm_cache_ != std::tie(velocity_cache, mass_cache, group_cache)) {
        LOG_TRACE("acquire centre-of-mass velocity");
        scoped_timer_type timer(runtime_.v_cm);
        std::tie(v_cm_, mean_mass_) = get_v_cm_and_mean_mass(*particle_, *group_);
        v_cm_cache_ = std::tie(velocity_cache, mass_cache, group_cache);
    }
    return v_cm_;
}

template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::vector_type const&
thermodynamics<dimension, float_type>::r_cm()
{
    cache<velocity_array_type> const& position_cache = particle_->position();
    cache<mass_array_type> const& mass_cache = particle_->mass();
    cache<size_type> const& group_cache = group_->size();

    if (r_cm_cache_ != std::tie(position_cache, mass_cache, group_cache)) {
        LOG_TRACE("acquire centre of mass");
        scoped_timer_type timer(runtime_.r_cm);
        r_cm_ = get_r_cm(*particle_, *group_, *box_);
        r_cm_cache_ = std::tie(position_cache, mass_cache, group_cache);
    }
    return r_cm_;
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::mean_mass()
{
    cache<velocity_array_type> const& velocity_cache = particle_->velocity();
    cache<mass_array_type> const& mass_cache = particle_->mass();
    cache<size_type> const& group_cache = group_->size();

    if (v_cm_cache_ != std::tie(velocity_cache, mass_cache, group_cache)) {
        LOG_TRACE("acquire mean particle mass");
        scoped_timer_type timer(runtime_.v_cm);
        std::tie(v_cm_, mean_mass_) = get_v_cm_and_mean_mass(*particle_, *group_);
        v_cm_cache_ = std::tie(velocity_cache, mass_cache, group_cache);
    }
    return mean_mass_;
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_pot()
{
    cache<en_pot_array_type> const& en_pot_cache = particle_->potential_energy();
    cache<size_type> const& group_cache = group_->size();

    if (en_pot_cache_ != std::tie(en_pot_cache, group_cache)) {
        LOG_TRACE("acquire potential energy");
        scoped_timer_type timer(runtime_.en_pot);
        en_pot_ = get_mean_en_pot(*particle_, *group_);
        en_pot_cache_ = std::tie(en_pot_cache, group_cache);
    }
    return en_pot_;
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::virial()
{
    cache<stress_pot_array_type> const& stress_pot_cache = particle_->stress_pot();
    cache<size_type> const& group_cache = group_->size();

    if (virial_cache_ != std::tie(stress_pot_cache, group_cache)) {
        LOG_TRACE("acquire virial");
        scoped_timer_type timer(runtime_.virial);
        virial_ = get_mean_virial(*particle_, *group_);
        virial_cache_ = std::tie(stress_pot_cache, group_cache);
    }
    return virial_;
}

template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::stress_tensor_type const&
thermodynamics<dimension, float_type>::stress_tensor()
{
    cache<stress_pot_array_type> const& stress_pot_cache = particle_->stress_pot();
    cache<velocity_array_type> const& velocity_cache = particle_->velocity();
    cache<size_type> const& group_cache = group_->size();

    if (stress_tensor_cache_ != std::tie(stress_pot_cache, velocity_cache, group_cache)) {
        LOG_TRACE("acquire stress tensor");
        scoped_timer_type timer(runtime_.stress_tensor);
        stress_tensor_ = get_stress_tensor(*particle_, *group_);
        stress_tensor_cache_ = std::tie(stress_pot_cache, velocity_cache, group_cache);
    }
    return stress_tensor_;
}

template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<thermodynamics, _Base>()
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("en_kin", &runtime::en_kin)
                        .def_readonly("force", &runtime::force)
                        .def_readonly("v_cm", &runtime::v_cm)
                        .def_readonly("r_cm", &runtime::r_cm)
                        .def_readonly("en_pot", &runtime::en_pot)
                        .def_readonly("virial", &runtime::virial)
                        .def_readonly("stress_tensor", &runtime::stress_tensor)
                ]
                .def_readonly("runtime", &thermodynamics::runtime_)

          , def("thermodynamics", &std::make_shared<thermodynamics
              , std::shared_ptr<particle_type>
              , std::shared_ptr<particle_group_type>
              , std::shared_ptr<box_type const>
              , std::shared_ptr<logger>
            >)
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

} // namespace host
} // namespace observables
} // namespace halmd
