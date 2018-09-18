/*
 * Copyright © 2015 Manuel Dibak
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

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <memory>

#include <halmd/mdsim/host/integrators/brownian.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>
#include <halmd/numeric/blas/detail/operators.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

// constructor
template <int dimension, typename float_type>
brownian<dimension, float_type>::brownian(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<random_type> random
  , std::shared_ptr<box_type const> box
  , double timestep
  , double T
  , matrix_type const& D
  , std::shared_ptr<logger> logger
)
  // dependency injection (initialize public variables)
  : particle_(particle)
  , random_(random)
  , box_(box)
  , D_(D)
  , logger_(logger)
{
    if (D_.size2() != 4) {
        throw std::invalid_argument("diffusion matrix has invalid shape: exactly 4 values per species are required");
    }

    set_timestep(timestep);
    set_temperature(T);

    LOG("diffusion constants: " << D_);
}

/**
 * set integration timestep
 */
template <int dimension, typename float_type>
void brownian<dimension, float_type>::set_timestep(double timestep)
{
    timestep_ = timestep;
    LOG("integration timestep: " << timestep_);
}

/**
 * set temperature of the heat bath
 */
template <int dimension, typename float_type>
void brownian<dimension, float_type>::set_temperature(double temperature)
{
    temperature_= temperature;
    LOG("temperature: " << temperature_);
}

/**
 * update both random and systematic parts of displacement
 */
template <int dimension, typename float_type>
void brownian<dimension, float_type>::update_displacement_(
    float_type const D
  , vector_type & r
  , vector_type & u
  , vector_type & v
  , vector_type const& f
)
{
    //random and systematic components of displacement
    vector_type dr_r, dr_s;

    //random part first
    float_type eta1, eta2;

    //not sure how sigma changes if D_perp != D_par
    float_type sigma = sqrt( 2 * timestep_ * D );
    std::tie(eta1, eta2) = random_->normal( sigma );
    dr_r[0] = eta1;
    dr_r[1] = eta2;
    if(dimension % 2) {
        float_type eta3, eta4;
        std::tie(eta3, eta4) = random_->normal( sigma );
        dr_r[2] = eta3;
    }

    //now systematic part (assumes D_perp = D_par)
    dr_s = v * timestep_ + D * timestep_ / temperature_ * f;

    r += dr_r + dr_s;

}

/**
 * free function to update orientation for 2d
 */
template<typename float_type>
void update_orientation_impl(
    float_type const D_rot
  , fixed_vector<float_type, 2> & u_2d
  , fixed_vector<float_type, 2> const& tau_2d
  , float_type timestep_
  , float_type temperature_
  , std::shared_ptr<random::host::random> random_
)
{
    //numerical limits for computation
    float_type epsilon = std::numeric_limits<float_type>::epsilon();

    //hack
    float_type eta1, eta2;
    fixed_vector<float_type, 3> u, tau, e1, e2;
    std::tie(u[0], u[1]) = std::tie(u_2d[0], u_2d[1]);
    std::tie(tau[0], tau[1]) = std::tie(tau_2d[0], tau_2d[1]);

    //construct trihedron along particle orientation
    if ( u[1] > epsilon || u[2] > epsilon) {
        e1[0] = 0; e1[1] = u[2]; e1[2] = -u[1];
    }
    else {
        e1[0] = u[1]; e1[1] = -u[0]; e1[2] = 0;
    }

    e1 /= norm_2(e1);
    e2  = cross_prod(u, e1);
    e2 /= norm_2(e2);
    float_type sigma_rot = sqrt( 2 * timestep_ * D_rot );
    std::tie(eta1, eta2) = random_->normal( sigma_rot );

    // first two terms are the random angular velocity, the final is the
    // systematic torque
    fixed_vector<float_type, 3> omega; 
    omega = eta1 * e1 + eta2 * e2 + tau * D_rot * timestep_ / temperature_ ;
    float_type  alpha        = norm_2(omega);
    omega              /= alpha;
    u = (1 - cos(alpha)) * inner_prod(omega, u) * omega + cos(alpha) * u + sin(alpha) * cross_prod(omega, u);
    u /= norm_2(u);

    // update the original 2d vector
    u_2d[0] = u[0];
    u_2d[1] = u[1];

}

/**
 * free function to update orientation for 3d
 */
template<typename float_type>
void update_orientation_impl(
    float_type const D_rot
  , fixed_vector<float_type, 3> & u
  , fixed_vector<float_type, 3> const& tau
  , float_type timestep_
  , float_type temperature_
  , std::shared_ptr<random::host::random> random_
)
{
    //numerical limits for computation
    float_type epsilon = std::numeric_limits<float_type>::epsilon();
    fixed_vector<float_type, 3>  e1, e2;
    float_type eta1, eta2;

    //construct trihedron along particle orientation
    if ( u[1] > epsilon || u[2] > epsilon) {
        e1[0] = 0; e1[1] = u[2]; e1[2] = -u[1];
    }
    else {
        e1[0] = u[1]; e1[1] = -u[0]; e1[2] = 0;
    }
    e2  = cross_prod(u, e1);

    // normalize vectors
    if ( (float) norm_2(e1) > 2e-38){
        e1 /= norm_2(e1);
    }
    if ( (float) norm_2(e2) > 2e-38){
        e2 /= norm_2(e2);
    }

    float_type sigma_rot = sqrt( 2 * timestep_ * D_rot );
    std::tie(eta1, eta2) = random_->normal( sigma_rot );

    fixed_vector<float_type, 3>  omega;

    // first two terms are the random angular velocity, the final is the
    // systematic torque
    omega = eta1 * e1 + eta2 * e2 + tau * D_rot * timestep_ / temperature_ ;
    float_type alpha = norm_2( omega );

    if ( (float) alpha > 2e-38){
        omega /= alpha;
    }

    u = (1 - cos(alpha)) * inner_prod(omega, u) * omega + cos(alpha) * u + sin(alpha) * cross_prod(omega, u);

    if ( (float) norm_2(u) > 2e-38){
        u /= norm_2(u);
    }
}

/**
 * Wrapper for the free functions that update the orientation
 */
template <int dimension, typename float_type>
void brownian<dimension, float_type>::update_orientation_(
    float_type const D_rot
  , vector_type & u
  , vector_type const& tau
)
{
   update_orientation_impl(D_rot, u, tau, timestep_, temperature_, random_);
}

/**
 * perform Brownian integration: update positions with random displacement 
 *
 * @f$ r(t + \Delta t) = \mu F(t) + \sigma d vec{W} @f$
 */
template <int dimension, typename float_type>
void brownian<dimension, float_type>::integrate()
{
    LOG_TRACE("update positions")

    //velocity_array_type const& velocity = read_cache(particle_->velocity());
    size_type nparticle = particle_->nparticle();
    auto const& torque  = read_cache(particle_->torque());
    auto const& force   = read_cache(particle_->force());


    // invalidate the particle caches after accessing the velocity!
    auto position       = make_cache_mutable(particle_->position());
    auto velocity       = make_cache_mutable(particle_->velocity());
    auto orientation    = make_cache_mutable(particle_->orientation());
    auto image          = make_cache_mutable(particle_->image());
    auto species        = make_cache_mutable(particle_->species());

    scoped_timer_type timer(runtime_.integrate);

    for (size_type i = 0 ; i < nparticle; ++i) {
        unsigned int particle_species = (*species)[i];
        float_type const D_perp    = D_(particle_species, 0);
        //float_type const D_par     = D_(particle_species, 1); 
        float_type const D_rot     = D_(particle_species, 2);
        //float_type prop_str  = D_(particle_species, 3);

        vector_type& r = (*position)[i];
        vector_type& v = (*velocity)[i];
        vector_type f  = force[i];
        vector_type& u = (*orientation)[i];
        vector_type tau = torque[i];

        update_displacement_(D_perp, r, u, v, f);

        // enforce periodic boundary conditions
        (*image)[i] += box_->reduce_periodic(r);
        
        // update orientation last (Ito interpretation)
        update_orientation_(D_rot, u, tau);
    }
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

template <int dimension, typename float_type>
void brownian<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<brownian>()
                    .property("integrate", &wrap_integrate<brownian>)
                    .property("timestep", &brownian::timestep)
                    .property("temperature", &brownian::temperature_)
                    .def("set_timestep", &brownian::set_timestep)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("integrate", &runtime::integrate)
                    ]
                    .def_readonly("runtime", &brownian::runtime_)

              , def("brownian", &std::make_shared<brownian
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<random_type>
                  , std::shared_ptr<box_type const>
                  , double
                  , double
                  , matrix_type const&
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_integrators_brownian(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    brownian<3, double>::luaopen(L);
    brownian<2, double>::luaopen(L);
#else
    brownian<3, float>::luaopen(L);
    brownian<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class brownian<3, double>;
template class brownian<2, double>;
#else
template class brownian<3, float>;
template class brownian<2, float>;
#endif

} // namespace integrators
} // namespace host
} // namespace mdsim
} // namespace halmd
