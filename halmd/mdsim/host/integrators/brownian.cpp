/*
 * Copyright © 2023 Jaslo Ziska
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

template <int dimension, typename float_type>
brownian<dimension, float_type>::brownian(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<random_type> random
  , std::shared_ptr<box_type const> box
  , float_type timestep
  , float_type temperature
  , matrix_type const& diff_const
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , random_(random)
  , box_(box)
  , diff_const_(diff_const)
  , logger_(logger)
{
    if (diff_const_.size1() != particle->nspecies()) {
        throw std::invalid_argument("diffusion matrix has invalid shape: exactly the number of species are required");
    }
    if (diff_const_.size2() != 2) {
        throw std::invalid_argument("diffusion matrix has invalid shape: exactly 2 values per species are required");
    }

    set_timestep(timestep);
    set_temperature(temperature);

    LOG("diffusion constants: " << diff_const_);
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
void brownian<dimension, float_type>::update_displacement(
    float_type const diff_const
  , vector_type& r
  , vector_type const& f
  , vector_type const& eta
)
{
    r += eta + diff_const * f * timestep_ / temperature_;
}

/**
 * free function to update orientation for 2d
 */
template<typename float_type>
void update_orientation_impl(
    float_type const diff_const
  , fixed_vector<float_type, 2>& u
  , fixed_vector<float_type, 2> const& tau
  , float_type timestep
  , float_type temperature
  , float_type eta
  , float_type ignore
)
{
    // numerical limits for computation
    float_type epsilon = std::numeric_limits<float_type>::epsilon();

    // the 2d case is just diffusion of orientation angle
    float_type omega = eta + tau[0] * diff_const * timestep / temperature;

    // rotate by omega
    fixed_vector<float_type, 2> e;
    e[0] = -u[1];
    e[1] = u[0];
    u = cos(omega) * u + sin(omega) * e;

    // ensure normalization, TODO: is this really necessary?
    if (norm_2(u) > epsilon) {
        u /= norm_2(u);
    }
}

/**
 * free function to update orientation for 3d
 */
template<typename float_type>
void update_orientation_impl(
    float_type const diff_const
  , fixed_vector<float_type, 3>& u
  , fixed_vector<float_type, 3> const& tau
  , float_type timestep
  , float_type temperature
  , float_type eta1
  , float_type eta2
)
{
    typedef fixed_vector<float_type, 3> vector_type;

    //numerical limits for computation
    float_type epsilon = std::numeric_limits<float_type>::epsilon();

    //construct trihedron along particle orientation
    vector_type e1, e2;
    if (u[1] > epsilon || u[2] > epsilon) {
        e1[0] = 0; e1[1] = u[2]; e1[2] = -u[1];
    } else {
        e1[0] = u[1]; e1[1] = -u[0]; e1[2] = 0;
    }
    e2  = cross_prod(u, e1);

    // normalize vectors, TODO: is this really necessary?
    if (norm_2(e1) > 2e-38){
        e1 /= norm_2(e1);
    }
    if (norm_2(e2) > 2e-38){
        e2 /= norm_2(e2);
    }

    // first two terms are the random angular velocity, the final is the
    // systematic torque
    vector_type omega = eta1 * e1 + eta2 * e2 + tau * diff_const * timestep / temperature;
    float_type alpha = norm_2(omega);

    if (alpha > 2e-38){
        omega /= alpha;
    }

    u = (1 - cos(alpha)) * inner_prod(omega, u) * omega + cos(alpha) * u + sin(alpha) * cross_prod(omega, u);

    if (norm_2(u) > 2e-38){
        u /= norm_2(u);
    }
}

/**
 * Wrapper for the free functions that update the orientation
 */
template <int dimension, typename float_type>
void brownian<dimension, float_type>::update_orientation(
    float_type const diff_const
  , vector_type& u
  , vector_type const& tau
  , float_type eta1
  , float_type eta2
)
{
    update_orientation_impl(diff_const, u, tau, timestep_, temperature_, eta1, eta2);
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

    size_type nparticle = particle_->nparticle();
    auto const& torque  = read_cache(particle_->torque());
    auto const& force   = read_cache(particle_->force());
    auto const& species = read_cache(particle_->species());

    // invalidate the particle caches after accessing the force and torque!
    auto position    = make_cache_mutable(particle_->position());
    auto orientation = make_cache_mutable(particle_->orientation());
    auto image       = make_cache_mutable(particle_->image());

    scoped_timer_type timer(runtime_.integrate);

    float_type rng_disp_cache = 0;
    bool rng_disp_cached = false;
    float_type rng_rot_cache = 0;
    bool rng_rot_cached = false;

    for (size_type i = 0 ; i < nparticle; ++i) {
        vector_type f = force[i];
        vector_type tau = torque[i];
        unsigned int s = species[i];

        vector_type& r = (*position)[i];
        vector_type& u = (*orientation)[i];

        float_type const diff_const_perp = diff_const_(s, 0);
        float_type const diff_const_rot = diff_const_(s, 1);
        float_type sigma_disp = sqrt(2 * timestep_ * diff_const_perp);
        float_type sigma_rot = sqrt(2 * timestep_ * diff_const_rot);

        // In the following do not generate the random numbers in the update_...() functions because we need to
        // cache the second random number (when an odd number of random numbers are required) for the next iteration of
        // the loop.

        // integrate positions
        vector_type dr;
        std::tie(dr[0], dr[1]) = random_->normal(sigma_disp);
        if (dimension % 2 == 1) {
            if (rng_disp_cached) {
                dr[2] = rng_disp_cache;
            } else {
                std::tie(dr[2], rng_disp_cache) = random_->normal(sigma_disp);
            }
            rng_disp_cached = !rng_disp_cached;
        }

        update_displacement(diff_const_perp, r, f, dr);

        // enforce periodic boundary conditions
        (*image)[i] += box_->reduce_periodic(r);

        // update orientation last (Ito interpretation)
        if (dimension == 2) {
            float_type eta;
            if (rng_rot_cached) {
                eta = rng_rot_cache;
            } else {
                std::tie(eta, rng_rot_cache) = random_->normal(sigma_rot);
            }
            rng_rot_cached = !rng_rot_cached;

            update_orientation(diff_const_rot, u, tau, eta, float_type(0));
        } else {
            float_type eta1, eta2;
            std::tie(eta1, eta2) = random_->normal(sigma_rot);

            update_orientation(diff_const_rot, u, tau, eta1, eta2);
        }
    }
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
                    .def("integrate", &brownian::integrate)
                    .def("set_timestep", &brownian::set_timestep)
                    .def("set_temperature", &brownian::set_temperature)
                    .property("timestep", &brownian::timestep)
                    .property("temperature", &brownian::temperature_)
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
