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
 * perform Brownian integration: update positions with random displacement
 *
 * @f$ r(t + \Delta t) = \mu F(t) + \sigma d vec{W} @f$
 */
template <int dimension, typename float_type>
void brownian<dimension, float_type>::integrate()
{
    LOG_TRACE("update positions")

    size_type nparticle = particle_->nparticle();
    auto const& force   = read_cache(particle_->force());
    auto const& species = read_cache(particle_->species());

    // invalidate the particle caches after accessing the force and torque!
    auto position    = make_cache_mutable(particle_->position());
    auto image       = make_cache_mutable(particle_->image());

    scoped_timer_type timer(runtime_.integrate);

    float_type rng_disp_cache = 0;
    bool rng_disp_cached = false;
    float_type rng_rot_cache = 0;
    bool rng_rot_cached = false;

    for (size_type i = 0 ; i < nparticle; ++i) {
        unsigned int s = species[i];

        // In the following do not generate the random numbers in the update_...() functions because we need to
        // cache the second random number (when an odd number of random numbers are required) for the next iteration of
        // the loop
        vector_type f = force[i];
        vector_type& r = (*position)[i];

        float_type const diff_const_perp = diff_const_(s, 0);
        float_type sigma_disp = sqrt(2 * timestep_ * diff_const_perp);

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
