/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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
#include <cmath>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/integrators/verlet_nvt_hoover.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using boost::fusion::at_key;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace integrators
{

template <int dimension, typename float_type>
verlet_nvt_hoover<dimension, float_type>::verlet_nvt_hoover(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , float_type timestep, float_type temperature
  , fixed_vector<float_type, 2>  const& mass
)
  // dependency injection
  : particle(particle)
  , box(box)
  // member initialisation
  , xi(0)
  , v_xi(0)
  , mass_xi_(mass)
{
    this->timestep(timestep);
    this->temperature(temperature);
    LOG("`masses' of heat bath variables: " << mass_xi_);
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
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

template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::temperature(double temperature)
{
    temperature_ = static_cast<float_type>(temperature);
    en_kin_target_2_ = dimension * particle->nbox * temperature_;

    LOG("temperature of heat bath: " << temperature_);
    LOG("target kinetic energy: " << en_kin_target_2_ / particle->nbox);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::integrate()
{
    scoped_timer<timer> timer_(at_key<integrate_>(runtime_));

    propagate_chain();

    for (size_t i = 0; i < particle->nbox; ++i) {
        vector_type& v = particle->v[i] += particle->f[i] * timestep_half_;
        vector_type& r = particle->r[i] += v * timestep_;
        // enforce periodic boundary conditions
        // TODO: reduction is now to (-L/2, L/2) instead of (0, L) as before
        // check that this is OK
        particle->image[i] += box->reduce_periodic(r);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::finalize()
{
    scoped_timer<timer> timer_(at_key<finalize_>(runtime_));

    // loop over all particles
    for (size_t i = 0; i < particle->nbox; ++i) {
        particle->v[i] += particle->f[i] * timestep_half_;
    }

    propagate_chain();
}

/**
 * propagate Nosé-Hoover chain
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::propagate_chain()
{
    scoped_timer<timer> timer_(at_key<propagate_>(runtime_));

    // compute total kinetic energy (multiplied by 2)
    float_type en_kin_2 = 0;
    BOOST_FOREACH(vector_type const& v, particle->v) {
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
    BOOST_FOREACH(vector_type& v, particle->v) {
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

template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L)
    [
        namespace_("halmd_wrapper")
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
                              , shared_ptr<box_type>
                              , float_type, float_type
                              , fixed_vector<float_type, 2> const&
                            >())
                            .def("register_runtimes", &verlet_nvt_hoover::register_runtimes)
                            .property("mass", &verlet_nvt_hoover::mass)
                            .property("module_name", &module_name_wrapper<dimension, float_type>)
                    ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &verlet_nvt_hoover<3, double>::luaopen
    ]
    [
        &verlet_nvt_hoover<2, double>::luaopen
    ];
#else
    [
        &verlet_nvt_hoover<3, float>::luaopen
    ]
    [
        &verlet_nvt_hoover<2, float>::luaopen
    ];
#endif
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
