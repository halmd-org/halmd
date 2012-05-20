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

#include <boost/make_shared.hpp>

#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/observables/gpu/thermodynamics.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
thermodynamics<dimension, float_type>::thermodynamics(
    boost::shared_ptr<particle_type const> particle
  , boost::shared_ptr<box_type const> box
  , boost::shared_ptr<clock_type const> clock
  , boost::shared_ptr<logger_type> logger
)
  // dependency injection
  : box_(box)
  , particle_(particle)
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
unsigned int thermodynamics<dimension, float_type>::nparticle() const
{
    return particle_->nparticle();
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
    if (!en_kin_.valid()) {
        LOG_TRACE("acquire kinetic energy");

        scoped_timer_type timer(runtime_.en_kin);
        en_kin_ = double(compute_en_kin_(particle_->velocity())()) / nparticle();
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
    if (!v_cm_.valid()) {
        LOG_TRACE("acquire centre-of-mass velocity");

        scoped_timer_type timer(runtime_.v_cm);
        v_cm_ = vector_type(compute_v_cm_(particle_->velocity())());
    }
    return v_cm_;
}

/**
 * compute potential energy
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_pot()
{
    if (!en_pot_.valid()) {
        LOG_TRACE("acquire potential energy");

        scoped_timer_type timer(runtime_.en_pot);
        en_pot_ = double(compute_en_pot_(particle_->en_pot())()) / nparticle();
    }
    return en_pot_;
}

/**
 * compute virial sum from potential part of stress tensor
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::virial()
{
    if (!virial_.valid()) {
        LOG_TRACE("acquire virial");

        scoped_timer_type timer(runtime_.virial);
        virial_ = double(compute_virial_(particle_->stress_pot())()) / nparticle();
    }
    return virial_;
}

/**
 * compute hypervirial sum
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::hypervirial()
{
    if (!hypervirial_.valid()) {
        LOG_TRACE("acquire hypervirial");

        scoped_timer_type timer(runtime_.hypervirial);
        hypervirial_ = double(compute_en_pot_(particle_->hypervirial())()) / nparticle();
    }
    return hypervirial_;
}

/**
 * clear data caches
 */
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
    static std::string class_name = demangled_name<thermodynamics>();
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                class_<thermodynamics, _Base>(class_name.c_str())
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
              , boost::shared_ptr<particle_type const>
              , boost::shared_ptr<box_type const>
              , boost::shared_ptr<clock_type const>
              , boost::shared_ptr<logger_type>
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

} // namespace observables
} // namespace gpu
} // namespace halmd
