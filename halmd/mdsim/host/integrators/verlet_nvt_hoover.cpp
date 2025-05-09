/*
 * Copyright © 2008-2014 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#include <halmd/config.hpp>

#include <algorithm>
#include <cmath>
#include <memory>

#include <halmd/mdsim/host/integrators/verlet_nvt_hoover.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

template <int dimension, typename float_type>
verlet_nvt_hoover<dimension, float_type>::verlet_nvt_hoover(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , float_type timestep
  , float_type temperature
  , float_type resonance_frequency
  , std::shared_ptr<logger> logger
)
  // public member initialisation
  : xi(0)
  , v_xi(0)
  // dependency injection
  , particle_(particle)
  , box_(box)
  , logger_(logger)
  // member initialisation
  , en_nhc_(0)
  , resonance_frequency_(resonance_frequency)
{
    set_timestep(timestep);

    LOG("resonance frequency of heat bath: " << resonance_frequency_);
    set_temperature(temperature);
}

template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::set_timestep(double timestep)
{
    timestep_ = timestep;
    timestep_half_ = timestep_ / 2;
    timestep_4_ = timestep_ / 4;
    timestep_8_ = timestep_ / 8;
}

/*
 * set temperature and adjust masses of heat bath variables
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::set_temperature(double temperature)
{
    temperature_ = static_cast<float_type>(temperature);
    en_kin_target_2_ = dimension * particle_->nparticle() * temperature_;

    // follow Martyna et al. [J. Chem. Phys. 97, 2635 (1992)]
    // for the masses of the heat bath variables
    float_type omega_sq = pow(2 * M_PI * resonance_frequency_, 2);
    unsigned int dof = dimension * particle_->nparticle();
    chain_type mass;
    mass[0] = dof * temperature_ / omega_sq;
    mass[1] = temperature_ / omega_sq;
    set_mass(mass);

    LOG("temperature of heat bath: " << temperature_);
    LOG_INFO("target kinetic energy: " << en_kin_target_2_ / particle_->nparticle());
}

template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::
set_mass(chain_type const& mass)
{
    mass_xi_ = mass;
    LOG_INFO("`mass' of heat bath variables: " << mass_xi_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::integrate()
{
    force_array_type const& force = read_cache(particle_->force());
    mass_array_type const& mass = read_cache(particle_->mass());
    size_type nparticle = particle_->nparticle();

    LOG_DEBUG("update positions and velocities: first leapfrog half-step")
    scoped_timer_type timer(runtime_.integrate);

    // invalidate the particle caches after accessing the force!
    auto position = make_cache_mutable(particle_->position());
    auto image = make_cache_mutable(particle_->image());
    auto velocity = make_cache_mutable(particle_->velocity());

    propagate_chain();

    for (size_type i = 0; i < nparticle; ++i) {
        vector_type& v = (*velocity)[i];
        vector_type& r = (*position)[i];
        v += force[i] * timestep_half_ / mass[i];
        r += v * timestep_;
        (*image)[i] += box_->reduce_periodic(r);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::finalize()
{
    force_array_type const& force = read_cache(particle_->force());
    mass_array_type const& mass = read_cache(particle_->mass());
    size_type nparticle = particle_->nparticle();

    LOG_DEBUG("update velocities: second leapfrog half-step")
    scoped_timer_type timer(runtime_.finalize);

    // invalidate the particle caches after accessing the force!
    auto velocity = make_cache_mutable(particle_->velocity());

    // loop over all particles
    for (size_type i = 0; i < nparticle; ++i) {
        (*velocity)[i] += force[i] * timestep_half_ / mass[i];
    }

    propagate_chain();

    // compute energy contribution of chain variables
    en_nhc_ = temperature_ * (dimension * nparticle * xi[0] + xi[1]);
    for (unsigned int i = 0; i < 2; ++i ) {
        en_nhc_ += mass_xi_[i] * v_xi[i] * v_xi[i] / 2;
    }
    en_nhc_ /= nparticle;
}

/**
 * propagate Nosé-Hoover chain
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::propagate_chain()
{
    auto velocity = make_cache_mutable(particle_->velocity());

    scoped_timer_type timer(runtime_.propagate);

    // compute total kinetic energy (multiplied by 2)
    float_type en_kin_2 = 0;
    for (vector_type const& v : *velocity) {
        // assuming unit mass for all particle types
        en_kin_2 += inner_prod(v, v);
    }

    // head of the chain
    v_xi[1] += (mass_xi_[0] * v_xi[0] * v_xi[0] - temperature_) / mass_xi_[1] * timestep_4_;
    float_type t = exp(-v_xi[1] * timestep_8_);
    v_xi[0] *= t;
    v_xi[0] += (en_kin_2 - en_kin_target_2_) / mass_xi_[0] * timestep_4_;
    v_xi[0] *= t;

    // propagate heat bath variables
    for (unsigned int i = 0; i < 2; ++i ) {
        xi[i] += v_xi[i] * timestep_half_;
    }

    // rescale velocities and kinetic energy
    float_type s = exp(-v_xi[0] * timestep_half_);
    for (vector_type& v : *velocity) {
        v *= s;
    }
    en_kin_2 *= s * s;

    // tail of the chain, mirrors the head
    v_xi[0] *= t;
    v_xi[0] += (en_kin_2 - en_kin_target_2_) / mass_xi_[0] * timestep_4_;
    v_xi[0] *= t;
    v_xi[1] += (mass_xi_[0] * v_xi[0] * v_xi[0] - temperature_) / mass_xi_[1] * timestep_4_;
}

template <typename integrator_type>
static std::function<void ()>
wrap_integrate(std::shared_ptr<integrator_type> self)
{
    return [=]() {
        self->integrate();
    };
}

template <typename integrator_type>
static std::function<void ()>
wrap_finalize(std::shared_ptr<integrator_type> self)
{
    return [=]() {
        self->finalize();
    };
}

template <typename integrator_type>
static std::function<typename integrator_type::chain_type& ()>
wrap_position(std::shared_ptr<integrator_type> self)
{
    return [=]() -> typename integrator_type::chain_type& {
        return self->xi;
    };
}

template <typename integrator_type>
static std::function<typename integrator_type::chain_type& ()>
wrap_velocity(std::shared_ptr<integrator_type> self)
{
    return [=]() -> typename integrator_type::chain_type& {
        return self->v_xi;
    };
}

template <typename integrator_type>
static std::function<double ()>
wrap_internal_energy(std::shared_ptr<integrator_type> self)
{
    return [=]() {
        return self->en_nhc();
    };
}

template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<verlet_nvt_hoover>()
                    .property("integrate", &wrap_integrate<verlet_nvt_hoover>)
                    .property("finalize", &wrap_finalize<verlet_nvt_hoover>)
                    .property("timestep", &verlet_nvt_hoover::timestep)
                    .property("temperature", &verlet_nvt_hoover::temperature)
                    .property("position", &wrap_position<verlet_nvt_hoover>)
                    .property("velocity", &wrap_velocity<verlet_nvt_hoover>)
                    .property("internal_energy", &wrap_internal_energy<verlet_nvt_hoover>)
                    .property("mass", &verlet_nvt_hoover::mass)
                    .property("resonance_frequency", &verlet_nvt_hoover::resonance_frequency)
                    .def("set_timestep", &verlet_nvt_hoover::set_timestep)
                    .def("set_temperature", &verlet_nvt_hoover::set_temperature)
                    .def("set_mass", &verlet_nvt_hoover::set_mass)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("integrate", &runtime::integrate)
                            .def_readonly("finalize", &runtime::finalize)
                            .def_readonly("propagate", &runtime::propagate)
                    ]
                    .def_readonly("runtime", &verlet_nvt_hoover::runtime_)

              , def("verlet_nvt_hoover", &std::make_shared<verlet_nvt_hoover
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type const>
                  , float_type
                  , float_type
                  , float_type
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_integrators_verlet_nvt_hoover(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    verlet_nvt_hoover<3, double>::luaopen(L);
    verlet_nvt_hoover<2, double>::luaopen(L);
#else
    verlet_nvt_hoover<3, float>::luaopen(L);
    verlet_nvt_hoover<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class verlet_nvt_hoover<3, double>;
template class verlet_nvt_hoover<2, double>;
#else
template class verlet_nvt_hoover<3, float>;
template class verlet_nvt_hoover<2, float>;
#endif

} // namespace integrators
} // namespace host
} // namespace mdsim
} // namespace halmd
