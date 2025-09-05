/*
 * Copyright © 2024 Max Orteu
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

#include <halmd/mdsim/host/integrators/brownian_euler.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>
#include <halmd/numeric/blas/detail/operators.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

template <int dimension, typename float_type>
brownian_euler<dimension, float_type>::brownian_euler(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<random_type> random
  , std::shared_ptr<box_type const> box
  , float_type timestep
  , float_type temperature
  , scalar_container_type const& diffusion
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , random_(random)
  , box_(box)
  , diffusion_(diffusion)
  , mobility_(diffusion.size())
  , noise_(diffusion.size())
  , logger_(logger)
{
    if (diffusion_.size() != particle_->nspecies()) {
        throw std::invalid_argument("diffusion constants have mismatching shape");
    }

    set_timestep(timestep);
    set_temperature(temperature);       // assigns mobility_

    for (size_t i = 0; i < diffusion_.size(); ++i) {
        noise_[i] = sqrt(2 * diffusion_[i]);
    }
    LOG_INFO("noise strengths: " << noise_);
    LOG("diffusion constants: " << diffusion_);
}

/**
 * set integration timestep
 */
template <int dimension, typename float_type>
void brownian_euler<dimension, float_type>::set_timestep(double timestep)
{
    timestep_ = timestep;
    LOG("integration timestep: " << timestep_);
}

/**
 * set temperature of the heat bath
 */
template <int dimension, typename float_type>
void brownian_euler<dimension, float_type>::set_temperature(double temperature)
{
    temperature_= temperature;
    LOG("temperature: " << temperature_);

    for (size_t i = 0; i < diffusion_.size(); ++i) {
        mobility_[i] = diffusion_[i] / temperature_;
    }
    LOG_INFO("mobility constants: " << mobility_);
}

/**
 * perform Brownian integration: update positions with random displacement
 *
 * @f$ r(t + \Delta t) = \mu F(t) + \sigma d vec{W} @f$
 */
template <int dimension, typename float_type>
void brownian_euler<dimension, float_type>::integrate()
{
    LOG_TRACE("update positions")

    size_type nparticle = particle_->nparticle();
    auto const& force   = read_cache(particle_->force());
    auto const& species = read_cache(particle_->species());

    // invalidate the particle caches only after accessing the force!
    auto position = make_cache_mutable(particle_->position());
    auto image    = make_cache_mutable(particle_->image());

    scoped_timer_type timer(runtime_.integrate);

    float_type rng_cache = 0;
    bool rng_cached = false;

    float_type sqrt_timestep_ = sqrt(timestep_);

    for (size_type i = 0 ; i < nparticle; ++i) {
        unsigned int s = species[i];

        // In the following do not generate the random numbers in the update_...() functions because we need to
        // cache the second random number (when an odd number of random numbers are required) for the next iteration of
        // the loop
        vector_type f = force[i];
        vector_type& r = (*position)[i];


        // draw Gaussian random vector
        vector_type dr;
        std::tie(dr[0], dr[1]) = random_->normal(float_type(1));
        if (dimension == 3) {
            if (rng_cached) {
                dr[2] = rng_cache;
            } else {
                std::tie(dr[2], rng_cache) = random_->normal(float_type(1));
            }
            rng_cached = !rng_cached;
        }
        dr *= noise_[s] * sqrt_timestep_;

        // integrate position: Euler-Maruyama scheme
        r += dr + (mobility_[s] * timestep_) * f;

        // enforce periodic boundary conditions
        (*image)[i] += box_->reduce_periodic(r);
    }
}

template <int dimension, typename float_type>
void brownian_euler<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<brownian_euler>()
                    .def("integrate", &brownian_euler::integrate)
                    .def("set_timestep", &brownian_euler::set_timestep)
                    .def("set_temperature", &brownian_euler::set_temperature)
                    .property("timestep", &brownian_euler::timestep)
                    .property("temperature", &brownian_euler::temperature_)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("integrate", &runtime::integrate)
                    ]
                    .def_readonly("runtime", &brownian_euler::runtime_)

              , def("brownian_euler", &std::make_shared<brownian_euler
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<random_type>
                  , std::shared_ptr<box_type const>
                  , double
                  , double
                  , scalar_container_type const&
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_integrators_brownian_euler(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    brownian_euler<3, double>::luaopen(L);
    brownian_euler<2, double>::luaopen(L);
#else
    brownian_euler<3, float>::luaopen(L);
    brownian_euler<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class brownian_euler<3, double>;
template class brownian_euler<2, double>;
#else
template class brownian_euler<3, float>;
template class brownian_euler<2, float>;
#endif

} // namespace integrators
} // namespace host
} // namespace mdsim
} // namespace halmd
