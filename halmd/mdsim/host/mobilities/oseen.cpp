/*
 * Copyright © 2011 Michael Kopp
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

#include <cmath> // sqrt(), pow(), M_PI
#include <algorithm> // fill

#include <halmd/mdsim/host/mobilities/oseen.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace mobilities {

template <int dimension, typename float_type>
oseen<dimension, float_type>::oseen(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , float radius
  , float viscosity
  , int order
  , boost::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
  // module parameters
  , radius_(radius)
  , viscosity_(viscosity)
  , order_(order)
{
    LOG("Particle radii: a = " << radius_ );
    LOG("Dynamic viscosity of fluid: η = " << viscosity_ );
    LOG("Order of accurancy of hydrodynamic interaction in (a/r): " << order_ );
    if( order_ <= 2 ) {
        LOG( "using Oseen tensor for hydrodynamic interaction" );
    }
    else if( order_ >= 3 ) {
        LOG( "using Rotne-Prager tensor for hydrodynamic interaction" );
    }
}

/**
 * \brief compute velocity from forces using Oseen Tensor calculus
 *
 * \note This algorithm exploits the fact that the Oseen Tensor is even in \f$
 * \vec r\f$ meaning that it computes to the same velocity regardless whether
 * \f$ \vec r\f$ or \f$ -\vec r\f$ is used.
 * This way \f$ r = \| \vec r \| \f$ needs only be computed \f$ N (N-1) \f$
 * times.
 *
 * \note The \a optimized code for interaction is taken from the GPU Module
 */
template <int dimension, typename float_type>
void oseen<dimension, float_type>::compute_velocity()
{
    scoped_timer<timer> timer_(runtime_.compute_velocity); // measure time 'till destruction

    // self mobility 1/(6 pi eta a)
    float_type self_mobility = 1 / (6 * M_PI * viscosity_ * radius_);

    for (unsigned int i = 0; i < particle_->nbox; ++i)
    {
        // self mobility
        particle_->v[i] += particle_->f[i];

        // interaction
        for (unsigned int j = i + 1; j < particle_->nbox; ++j)
        {
            // vector connecting the two particles i and j
            vector_type dr =  particle_->r[i] - particle_->r[j];
            // apply minimum image convention in PBC
            box_->reduce_periodic( dr );
            // distance between particles
            float_type dist2 = inner_prod(dr, dr);
            float_type dist = sqrt( dist2 );
            float_type b = radius_ / dist; //< to simplify following calculations

            if (order_ <= 2) { //oseen
                particle_->v[i] += (particle_->f[j] + (inner_prod(dr,particle_->f[j]) / dist2) * dr) * 0.75f * b;
                particle_->v[j] += (particle_->f[i] + (inner_prod(dr,particle_->f[i]) / dist2) * dr) * 0.75f * b;
            }
            else if (order_ <= 4) { // rotne prager
                if (dist < 2 * radius_) { // close branch
                    LOG_ONCE( "Particles are at distance " << dist << " -- using close branch" );
                    particle_->v[i] += ( 1 - (9.f / 32) * dist / radius_ ) * particle_->f[j] + ( (3.f / 32) * inner_prod(dr, particle_->f[j]) / (radius_ * dist) ) * dr;
                    particle_->v[j] += ( 1 - (9.f / 32) * dist / radius_ ) * particle_->f[i] + ( (3.f / 32) * inner_prod(dr, particle_->f[i]) / (radius_ * dist) ) * dr;
                }
                else { // default branch
                    float_type b2 = b * b;
                    particle_->v[i] += ((0.75f + 0.5f * b2) * b) * particle_->f[j] + ((0.75f - 1.5f * b2) * b * inner_prod(dr, particle_->f[j]) / dist2) * dr;
                    particle_->v[j] += ((0.75f + 0.5f * b2) * b) * particle_->f[i] + ((0.75f - 1.5f * b2) * b * inner_prod(dr, particle_->f[i]) / dist2) * dr;
                }
            }
        }
        particle_->v[i] *= self_mobility; //< this has been factorized in previous computations
    }
}

//! compute oseen tensor -- NOT YET IMPLEMENTED
template <int dimension, typename float_type>
void oseen<dimension, float_type>::compute() {
    scoped_timer<timer> timer_(runtime_.compute); // measure time 'till destruction
}

// Wrappers expose signal-functions which can passed to a signal.
template <typename mobility_type>
typename signal<void ()>::slot_function_type
wrap_compute(shared_ptr<mobility_type> self)
{
    return bind(&mobility_type::compute, self);
}

template <typename mobility_type>
typename signal<void ()>::slot_function_type
wrap_compute_velocity(shared_ptr<mobility_type> self)
{
    return bind(&mobility_type::compute_velocity, self);
}

//! register class in luabind
template <int dimension, typename float_type>
void oseen<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("oseen_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("mobilities")
                [
                    class_<oseen, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type>
                          , shared_ptr<box_type>
                          , double
                          , double
                          , int
                          , shared_ptr<logger>
                        >())
                        .property("radius", &oseen::radius)
                        .property("viscosity", &oseen::viscosity)
                        .property("order", &oseen::order)
                        .property("compute", &wrap_compute<oseen>)
                        .property("compute_velocity", &wrap_compute_velocity<oseen>)
                        // register runtime accumulators
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("compute", &runtime::compute)
                                .def_readonly("compute_velocity", &runtime::compute_velocity)
                        ]
                        .def_readonly("runtime", &oseen::runtime_)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_mobilities_oseen(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    oseen<3, double>::luaopen(L);
    oseen<2, double>::luaopen(L);
#else
    oseen<3, float>::luaopen(L);
    oseen<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class oseen<3, double>;
template class oseen<2, double>;
#else
template class oseen<3, float>;
template class oseen<2, float>;
#endif

} // namespace mobilities
} // namespace host
} // namespace mdsim
} // namespace halmd
