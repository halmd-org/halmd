/*
 * Copyright © 2010-2012 Felix Höfling
 * Copyright © 2010-2012 Peter Colberg
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

#include <memory>

#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/observables/gpu/thermodynamics.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
thermodynamics<dimension, float_type>::thermodynamics(
    std::shared_ptr<particle_type const> particle
  , std::shared_ptr<force_type> force
  , std::shared_ptr<particle_group_type> group
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<logger_type> logger
)
  : particle_(particle)
  , force_(force)
  , group_(group)
  , box_(box)
  , logger_(logger)
{
}

template <int dimension, typename float_type>
unsigned int thermodynamics<dimension, float_type>::nparticle() const
{
    return read_cache(group_->size());
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::volume() const
{
    return box_->volume();
}

/**
 * compute mean kinetic energy per particle
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_kin()
{
    cache<velocity_array_type> const& velocity_cache = particle_->velocity();
    cache<size_type> const& group_cache = group_->size();

    if (en_kin_cache_ != std::tie(velocity_cache, group_cache)) {
        LOG_TRACE("acquire kinetic energy");
        scoped_timer_type timer(runtime_.en_kin);
        en_kin_ = get_mean_en_kin(*particle_, *group_);
        en_kin_cache_ = std::tie(velocity_cache, group_cache);
    }
    return en_kin_;
}

/**
 * compute mean velocity
 */
template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::vector_type const&
thermodynamics<dimension, float_type>::v_cm()
{
    cache<velocity_array_type> const& velocity_cache = particle_->velocity();
    cache<size_type> const& group_cache = group_->size();

    if (v_cm_cache_ != std::tie(velocity_cache, group_cache)) {
        LOG_TRACE("acquire centre-of-mass velocity");
        scoped_timer_type timer(runtime_.v_cm);
        std::tie(v_cm_, mean_mass_) = get_v_cm_and_mean_mass(*particle_, *group_);
        v_cm_cache_ = std::tie(velocity_cache, group_cache);
    }
    return v_cm_;
}

/**
 * compute centre of mass
 */
template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::vector_type const&
thermodynamics<dimension, float_type>::r_cm()
{
    cache<velocity_array_type> const& position_cache = particle_->position();
    cache<size_type> const& group_cache = group_->size();

    if (r_cm_cache_ != std::tie(position_cache, group_cache)) {
        LOG_TRACE("acquire centre-of-mass velocity");
        scoped_timer_type timer(runtime_.v_cm);
        r_cm_ = get_r_cm(*particle_, *group_, *box_);
        r_cm_cache_ = std::tie(position_cache, group_cache);
    }
    return r_cm_;
}

/**
 * compute mean particle mass
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::mean_mass()
{
    cache<velocity_array_type> const& velocity_cache = particle_->velocity();
    cache<size_type> const& group_cache = group_->size();

    if (v_cm_cache_ != std::tie(velocity_cache, group_cache)) {
        LOG_TRACE("acquire centre-of-mass velocity");
        scoped_timer_type timer(runtime_.v_cm);
        std::tie(v_cm_, mean_mass_) = get_v_cm_and_mean_mass(*particle_, *group_);
        v_cm_cache_ = std::tie(velocity_cache, group_cache);
    }
    return mean_mass_;
}

/**
 * compute potential energy
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_pot()
{
    cache<en_pot_array_type> const& en_pot_cache = force_->en_pot();
    cache<size_type> const& group_cache = group_->size();

    if (en_pot_cache_ != std::tie(en_pot_cache, group_cache)) {
        LOG_TRACE("acquire potential energy");
        scoped_timer_type timer(runtime_.en_pot);
        en_pot_ = get_mean_en_pot(*force_, *group_);
        en_pot_cache_ = std::tie(en_pot_cache, group_cache);
    }
    return en_pot_;
}

/**
 * compute virial sum from potential part of stress tensor
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::virial()
{
    cache<stress_pot_array_type> const& stress_pot_cache = force_->stress_pot();
    cache<size_type> const& group_cache = group_->size();

    if (virial_cache_ != std::tie(stress_pot_cache, group_cache)) {
        LOG_TRACE("acquire virial");
        scoped_timer_type timer(runtime_.virial);
        virial_ = get_mean_virial(*force_, *group_);
        virial_cache_ = std::tie(stress_pot_cache, group_cache);
    }
    return virial_;
}

/**
 * compute hypervirial sum
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::hypervirial()
{
    cache<hypervirial_array_type> const& hypervirial_cache = force_->hypervirial();
    cache<size_type> const& group_cache = group_->size();

    if (hypervirial_cache_ != std::tie(hypervirial_cache, group_cache)) {
        LOG_TRACE("acquire hypervirial");
        scoped_timer_type timer(runtime_.hypervirial);
        hypervirial_ = get_mean_hypervirial(*force_, *group_);
        hypervirial_cache_ = std::tie(hypervirial_cache, group_cache);
    }
    return hypervirial_;
}

template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                class_<thermodynamics, _Base>()
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("en_kin", &runtime::en_kin)
                            .def_readonly("v_cm", &runtime::v_cm)
                            .def_readonly("r_cm", &runtime::r_cm)
                            .def_readonly("en_pot", &runtime::en_pot)
                            .def_readonly("virial", &runtime::virial)
                            .def_readonly("hypervirial", &runtime::hypervirial)
                    ]
                    .def_readonly("runtime", &thermodynamics::runtime_)
            ]
        ]

      , namespace_("observables")
        [
            def("thermodynamics", &std::make_shared<thermodynamics
              , std::shared_ptr<particle_type const>
              , std::shared_ptr<force_type>
              , std::shared_ptr<particle_group_type>
              , std::shared_ptr<box_type const>
              , std::shared_ptr<logger_type>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_thermodynamics(lua_State* L)
{
    thermodynamics<3, float>::luaopen(L);
    thermodynamics<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class thermodynamics<3, float>;
template class thermodynamics<2, float>;

} // namespace gpu
} // namespace observables
} // namespace halmd
