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
#include <boost/bind.hpp>
#include <cmath>
#include <memory>

#include <halmd/mdsim/host/integrators/verlet_nvt_andersen.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

template <int dimension, typename float_type>
verlet_nvt_andersen<dimension, float_type>::verlet_nvt_andersen(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<random_type> random
  , float_type timestep
  , float_type temperature
  , float_type coll_rate
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , box_(box)
  , random_(random)
  , coll_rate_(coll_rate)
  , logger_(logger)
{
    set_timestep(timestep);
    set_temperature(temperature);
    LOG("collision rate with heat bath: " << coll_rate_);
}

template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::set_timestep(double timestep)
{
    timestep_ = timestep;
    timestep_half_ = 0.5 * timestep;
    coll_prob_ = coll_rate_ * timestep;
}

template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::set_temperature(double temperature)
{
    temperature_ = temperature;
    sqrt_temperature_ = std::sqrt(temperature_);
    LOG("temperature of heat bath: " << temperature_);
}

template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::integrate()
{
    LOG_TRACE("update positions and velocities")

    force_array_type const& force = read_cache(particle_->force());
    mass_array_type const& mass = read_cache(particle_->mass());
    size_type nparticle = particle_->nparticle();

    // invalidate the particle caches after accessing the force!
    auto position = make_cache_mutable(particle_->position());
    auto image = make_cache_mutable(particle_->image());
    auto velocity = make_cache_mutable(particle_->velocity());

    scoped_timer_type timer(runtime_.integrate);

    for (size_type i = 0; i < nparticle; ++i) {
        vector_type& v = (*velocity)[i];
        vector_type& r = (*position)[i];
        v += force[i] * timestep_half_ / mass[i];
        r += v * timestep_;
        (*image)[i] += box_->reduce_periodic(r);
    }
}

template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::finalize()
{
    LOG_TRACE("update velocities")

    force_array_type const& force = read_cache(particle_->force());
    mass_array_type const& mass = read_cache(particle_->mass());
    size_type nparticle = particle_->nparticle();

    // invalidate the particle caches after accessing the force!
    auto velocity = make_cache_mutable(particle_->velocity());

    scoped_timer_type timer(runtime_.finalize);

    // cache random numbers
    float_type rng_cache = 0;
    bool rng_cache_valid = false;

    // loop over all particles
    for (size_type i = 0; i < nparticle; ++i) {
        vector_type& v = (*velocity)[i];
        // is deterministic step?
        if (random_->uniform<float_type>() > coll_prob_) {
            v += force[i] * timestep_half_ / mass[i];
        }
        // stochastic coupling with heat bath
        else {
            // assign two velocity components at a time
            for (unsigned int i = 0; i < dimension - 1; i += 2) {
                std::tie(v[i], v[i + 1]) = random_->normal(sqrt_temperature_);
            }
            // handle last component separately for odd dimensions
            if (dimension % 2 == 1) {
                if (rng_cache_valid) {
                    v[dimension - 1] = rng_cache;
                }
                else {
                    std::tie(v[dimension - 1], rng_cache) = random_->normal(sqrt_temperature_);
                }
                rng_cache_valid = !rng_cache_valid;
            }
        }
    }
}

template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<verlet_nvt_andersen>()
                    .def("integrate", &verlet_nvt_andersen::integrate)
                    .def("finalize", &verlet_nvt_andersen::finalize)
                    .def("set_timestep", &verlet_nvt_andersen::set_timestep)
                    .def("set_temperature", &verlet_nvt_andersen::set_temperature)
                    .property("timestep", &verlet_nvt_andersen::timestep)
                    .property("temperature", &verlet_nvt_andersen::temperature)
                    .property("collision_rate", &verlet_nvt_andersen::collision_rate)
                    .scope
                    [
                        class_<runtime>()
                            .def_readonly("integrate", &runtime::integrate)
                            .def_readonly("finalize", &runtime::finalize)
                    ]
                    .def_readonly("runtime", &verlet_nvt_andersen::runtime_)

              , def("verlet_nvt_andersen", &std::make_shared<verlet_nvt_andersen
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type const>
                  , std::shared_ptr<random_type>
                  , float_type
                  , float_type
                  , float_type
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_integrators_verlet_nvt_andersen(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    verlet_nvt_andersen<3, double>::luaopen(L);
    verlet_nvt_andersen<2, double>::luaopen(L);
#else
    verlet_nvt_andersen<3, float>::luaopen(L);
    verlet_nvt_andersen<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class verlet_nvt_andersen<3, double>;
template class verlet_nvt_andersen<2, double>;
#else
template class verlet_nvt_andersen<3, float>;
template class verlet_nvt_andersen<2, float>;
#endif

} // namespace integrators
} // namespace host
} // namespace mdsim
} // namespace halmd
