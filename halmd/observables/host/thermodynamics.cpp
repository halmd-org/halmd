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

#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
thermodynamics<dimension, float_type>::thermodynamics(
    std::shared_ptr<particle_type const> particle
  , std::shared_ptr<particle_group_type> group
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , group_(group)
  , box_(box)
  , logger_(logger)
{
}

template <int dimension, typename float_type>
unsigned int thermodynamics<dimension, float_type>::nparticle() const
{
    cache_proxy<size_type const> size = group_->size();
    return *size;
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

        cache_proxy<group_array_type const> group = group_->unordered();
        cache_proxy<velocity_array_type const> velocity = velocity_cache;
        cache_proxy<mass_array_type const> mass = mass_cache;

        scoped_timer_type timer(runtime_.en_kin);

        double mv2 = 0;
        for (size_type i : *group) {
            mv2 += (*mass)[i] * inner_prod((*velocity)[i], (*velocity)[i]);
        }
        en_kin_ = 0.5 * mv2 / group->size();

        en_kin_cache_ = std::tie(velocity_cache, mass_cache, group_cache);
    }
    return en_kin_;
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

        cache_proxy<group_array_type const> group = group_->unordered();
        cache_proxy<velocity_array_type const> velocity = velocity_cache;
        cache_proxy<mass_array_type const> mass = mass_cache;

        scoped_timer_type timer(runtime_.v_cm);

        vector_type mv = 0;
        double m = 0;
        for (size_type i : *group) {
            mv += (*mass)[i] * (*velocity)[i];
            m += (*mass)[i];
        }
        v_cm_ = mv / m;

        v_cm_cache_ = std::tie(velocity_cache, mass_cache, group_cache);
    }
    return v_cm_;
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_pot()
{
    cache<en_pot_array_type> const& en_pot_cache = particle_->en_pot();
    cache<size_type> const& group_cache = group_->size();

    if (en_pot_cache_ != std::tie(en_pot_cache, group_cache)) {
        LOG_TRACE("acquire potential energy");

        cache_proxy<group_array_type const> group = group_->unordered();
        cache_proxy<en_pot_array_type const> en_pot = en_pot_cache;

        scoped_timer_type timer(runtime_.en_pot);

        double sum = 0;
        for (size_type i : *group) {
            sum += (*en_pot)[i];
        }
        en_pot_ = sum / group->size();

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

        cache_proxy<group_array_type const> group = group_->unordered();
        cache_proxy<stress_pot_array_type const> stress_pot = stress_pot_cache;

        scoped_timer_type timer(runtime_.virial);

        double sum = 0;
        for (size_type i : *group) {
            sum += (*stress_pot)[i][0];
        }
        virial_ = sum / group->size();

        virial_cache_ = std::tie(stress_pot_cache, group_cache);
    }
    return virial_;
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::hypervirial()
{
    cache<hypervirial_array_type> const& hypervirial_cache = particle_->hypervirial();
    cache<size_type> const& group_cache = group_->size();

    if (hypervirial_cache_ != std::tie(hypervirial_cache, group_cache)) {
        LOG_TRACE("acquire hypervirial");

        cache_proxy<group_array_type const> group = group_->unordered();
        cache_proxy<hypervirial_array_type const> hypervirial = hypervirial_cache;

        scoped_timer_type timer(runtime_.hypervirial);

        double sum = 0;
        for (size_type i : *group) {
            sum += (*hypervirial)[i];
        }
        hypervirial_ = sum / group->size();

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
            class_<thermodynamics, _Base>()
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("en_kin", &runtime::en_kin)
                        .def_readonly("v_cm", &runtime::v_cm)
                        .def_readonly("en_pot", &runtime::en_pot)
                        .def_readonly("virial", &runtime::virial)
                        .def_readonly("hypervirial", &runtime::hypervirial)
                ]
                .def_readonly("runtime", &thermodynamics::runtime_)

          , def("thermodynamics", &std::make_shared<thermodynamics
              , std::shared_ptr<particle_type const>
              , std::shared_ptr<particle_group_type>
              , std::shared_ptr<box_type const>
              , std::shared_ptr<logger_type>
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
