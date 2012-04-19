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

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <cmath>
#include <string>

#include <halmd/mdsim/host/integrators/verlet_nvt_hoover.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace integrators
{

template <int dimension, typename float_type>
verlet_nvt_hoover<dimension, float_type>::verlet_nvt_hoover(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type const> box
  , float_type timestep
  , float_type temperature
  , float_type resonance_frequency
  , shared_ptr<logger_type> logger
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
    this->timestep(timestep);

    LOG("resonance frequency of heat bath: " << resonance_frequency_);
    this->temperature(temperature);
}

template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::timestep(double timestep)
{
    timestep_ = static_cast<float_type>(timestep);
    timestep_half_ = timestep_ / 2;
    timestep_4_ = timestep_ / 4;
    timestep_8_ = timestep_ / 8;

    LOG("integration timestep: " << timestep_);
}

/*
 * set temperature and adjust masses of heat bath variables
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::temperature(double temperature)
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
    LOG_DEBUG("target kinetic energy: " << en_kin_target_2_ / particle_->nparticle());
}

template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::
set_mass(chain_type const& mass)
{
    mass_xi_ = mass;
    LOG("`mass' of heat bath variables: " << mass_xi_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::integrate()
{
    scoped_timer_type timer(runtime_.integrate);

    typename particle_type::force_array_type& force = particle_->force();

    propagate_chain();

    for (size_t i = 0; i < particle_->nparticle(); ++i) {
        vector_type& v = particle_->v[i] += force[i] * timestep_half_;
        vector_type& r = particle_->r[i] += v * timestep_;
        // enforce periodic boundary conditions
        particle_->image[i] += box_->reduce_periodic(r);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::finalize()
{
    scoped_timer_type timer(runtime_.finalize);

    typename particle_type::force_array_type& force = particle_->force();

    // loop over all particles
    for (size_t i = 0; i < particle_->nparticle(); ++i) {
        particle_->v[i] += force[i] * timestep_half_;
    }

    propagate_chain();

    // compute energy contribution of chain variables
    en_nhc_ = temperature_ * (dimension * particle_->nparticle() * xi[0] + xi[1]);
    for (unsigned int i = 0; i < 2; ++i ) {
        en_nhc_ += mass_xi_[i] * v_xi[i] * v_xi[i] / 2;
    }
    en_nhc_ /= particle_->nparticle();
}

/**
 * propagate Nosé-Hoover chain
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::propagate_chain()
{
    scoped_timer_type timer(runtime_.propagate);

    // compute total kinetic energy (multiplied by 2)
    float_type en_kin_2 = 0;
    BOOST_FOREACH(vector_type const& v, particle_->v) {
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
    BOOST_FOREACH(vector_type& v, particle_->v) {
        v *= s;
    }
    en_kin_2 *= s * s;

    // tail of the chain, mirrors the head
    v_xi[0] *= t;
    v_xi[0] += (en_kin_2 - en_kin_target_2_) / mass_xi_[0] * timestep_4_;
    v_xi[0] *= t;
    v_xi[1] += (mass_xi_[0] * v_xi[0] * v_xi[0] - temperature_) / mass_xi_[1] * timestep_4_;
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(verlet_nvt_hoover<dimension, float_type> const&)
{
    return verlet_nvt_hoover<dimension, float_type>::module_name();
}

template <typename integrator_type>
static function<typename integrator_type::chain_type& ()>
wrap_xi(shared_ptr<integrator_type> integrator)
{
    return lambda::var(integrator->xi);
}

template <typename integrator_type>
static function<typename integrator_type::chain_type& ()>
wrap_v_xi(shared_ptr<integrator_type> integrator)
{
    return lambda::var(integrator->v_xi);
}

template <typename integrator_type>
static function<double ()>
wrap_en_nhc(shared_ptr<integrator_type> integrator)
{
    return bind(&integrator_type::en_nhc, integrator);
}

template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("integrators")
                [
                    class_<
                        verlet_nvt_hoover
                      , shared_ptr<_Base_Base>
                      , bases<_Base_Base, _Base>
                    >(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type>
                          , shared_ptr<box_type const>
                          , float_type, float_type, float_type
                          , shared_ptr<logger_type>
                        >())
                        .property("xi", &wrap_xi<verlet_nvt_hoover>)
                        .property("v_xi", &wrap_v_xi<verlet_nvt_hoover>)
                        .property("en_nhc", &wrap_en_nhc<verlet_nvt_hoover>)
                        .property("mass", &verlet_nvt_hoover::mass)
                        .property("resonance_frequency", &verlet_nvt_hoover::resonance_frequency)
                        .property("module_name", &module_name_wrapper<dimension, float_type>)
                        .def("set_mass", &verlet_nvt_hoover::set_mass)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("integrate", &runtime::integrate)
                                .def_readonly("finalize", &runtime::finalize)
                                .def_readonly("propagate", &runtime::propagate)
                        ]
                        .def_readonly("runtime", &verlet_nvt_hoover::runtime_)
                ]
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

}}} // namespace mdsim::host::integrators

} // namespace halmd
