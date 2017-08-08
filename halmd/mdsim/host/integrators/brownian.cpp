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
#include <cmath>
#include <memory>

#include <halmd/mdsim/host/integrators/brownian.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

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
  , temperature_(T)
  , D_(D)
  , logger_(logger)
{
    set_timestep(timestep);
}

/**
 * set integration timestep
 */
template <int dimension, typename float_type>
void brownian<dimension, float_type>::set_timestep(double timestep)
{
    timestep_ = timestep;
}

/**
 * set temperature of the heat bath
 */
template <int dimension, typename float_type>
void brownian<dimension, float_type>::set_temperature(double temperature)
{
    temperature_= temperature;
}

/**
 * compute random displacement
 */

template <int dimension, typename float_type>
typename brownian<dimension, float_type>::vector_type brownian<dimension, float_type>::random_displacement_(double D)
{
   float_type eta1, eta2;
   vector_type dr;
   float_type sigma = sqrt( 2 * timestep_ * D );
   std::tie(eta1, eta2) = random_->normal( sigma );
   dr[0] = eta1;
   dr[1] = eta2;
   if(dimension % 2) {
       float_type eta3, eta4;
       std::tie(eta3, eta4) = random_->normal( sigma );
       dr[2] = eta3;
   }
   return dr;
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

    // invalidate the particle caches after accessing the velocity!
    auto position = make_cache_mutable(particle_->position());
    auto image = make_cache_mutable(particle_->image());
    auto species = make_cache_mutable(particle_->species());
    scoped_timer_type timer(runtime_.integrate);
    float_type D;
    unsigned int particle_species;

    for (size_type i = 0 ; i < nparticle; ++i) {
        particle_species = (*species)[i];
        D = D_.data()[ particle_species ];
        vector_type& r = (*position)[i];
        vector_type dr = random_displacement_(D);
        r += dr;
        (*image)[i] += box_->reduce_periodic(r);
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
