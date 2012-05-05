/*
 * Copyright © 2010-2012  Felix Höfling and Peter Colberg
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

#include <boost/foreach.hpp>

#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
thermodynamics<dimension, float_type>::thermodynamics(
    boost::shared_ptr<particle_group_type const> particle_group
  , boost::shared_ptr<box_type const> box
  , boost::shared_ptr<clock_type const> clock
  , boost::shared_ptr<logger_type> logger
)
  // dependency injection
  : box_(box)
  , particle_group_(particle_group)
  , particle_(particle_group->particle())
  , logger_(logger)
  // initialise members
  , en_kin_(clock)
  , v_cm_(clock)
  , en_pot_(clock)
  , virial_(clock)
  , hypervirial_(clock)
{
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_kin()
{
    if (!en_kin_.valid()) {
        LOG_TRACE("acquire kinetic energy");

        scoped_timer_type timer(runtime_.en_kin);

        typename particle_type::velocity_array_type const& velocity = particle_->velocity();
        typename particle_type::mass_array_type const& mass = particle_->mass();

        // FIXME use particle_group_->selection_mask()
        double mv2 = 0;
        for (size_t i = 0; i < particle_->nparticle(); ++i) {
            // assuming unit mass for all particle types
            mv2 += mass[i] * inner_prod(velocity[i], velocity[i]);
        }
        en_kin_ = 0.5 * mv2 / nparticle();
    }
    return en_kin_;
}

template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::vector_type const&
thermodynamics<dimension, float_type>::v_cm()
{
    if (!v_cm_.valid()) {
        LOG_TRACE("acquire centre-of-mass velocity");

        scoped_timer_type timer(runtime_.v_cm);

        typename particle_type::velocity_array_type const& velocity = particle_->velocity();
        typename particle_type::mass_array_type const& mass = particle_->mass();

        // FIXME use particle_group_->selection_mask()
        vector_type mv = 0;
        double m = 0;
        for (size_t i = 0; i < particle_->nparticle(); ++i) {
            mv += mass[i] * velocity[i];
            m += mass[i];
        }
        v_cm_ = mv / m;
    }
    return v_cm_;
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_pot()
{
    typedef typename particle_type::en_pot_type en_pot_type;
    typename particle_type::en_pot_array_type const& en_pot = particle_->en_pot();

    if (!en_pot_.valid()) {
        LOG_TRACE("acquire potential energy");

        scoped_timer_type timer(runtime_.en_pot);

        // FIXME use particle_group_->selection_mask()
        double sum = 0;
        BOOST_FOREACH(en_pot_type value, en_pot) {
            sum += value;
        }
        en_pot_ = sum / nparticle();
    }
    return en_pot_;
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::virial()
{
    typedef typename particle_type::stress_pot_type stress_pot_type;
    typename particle_type::stress_pot_array_type const& stress_pot = particle_->stress_pot();

    if (!virial_.valid()) {
        LOG_TRACE("acquire virial");

        scoped_timer_type timer(runtime_.virial);

        // FIXME use particle_group_->selection_mask()
        double sum = 0;
        BOOST_FOREACH(stress_pot_type const& value, stress_pot) {
            sum += value[0];
        }
        virial_ = sum / nparticle();
    }
    return virial_;
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::hypervirial()
{
    typedef typename particle_type::hypervirial_type hypervirial_type;
    typename particle_type::hypervirial_array_type const& hypervirial = particle_->hypervirial();

    if (!hypervirial_.valid()) {

        scoped_timer_type timer(runtime_.hypervirial);

        LOG_TRACE("acquire hypervirial");
        // FIXME use particle_group_->selection_mask()
        double sum = 0;
        BOOST_FOREACH(hypervirial_type value, hypervirial) {
            sum += value;
        }
        hypervirial_ = sum / nparticle();
    }
    return hypervirial_;
}

template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::clear_cache()
{
    en_kin_.clear();
    v_cm_.clear();
    en_pot_.clear();
    virial_.clear();
    hypervirial_.clear();
}

template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("thermodynamics_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<thermodynamics, boost::shared_ptr<_Base>, _Base>(class_name.c_str())
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
            ]
        ]

      , namespace_("observables")
        [
            def("thermodynamics", &boost::make_shared<thermodynamics
              , boost::shared_ptr<particle_group_type const>
              , boost::shared_ptr<box_type const>
              , boost::shared_ptr<clock_type const>
              , boost::shared_ptr<logger_type>
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

} // namespace observables
} // namespace host
} // namespace halmd
