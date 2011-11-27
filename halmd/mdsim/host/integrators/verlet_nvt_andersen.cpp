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
#include <cmath>
#include <string>

#include <halmd/mdsim/host/integrators/verlet_nvt_andersen.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

template <int dimension, typename float_type>
verlet_nvt_andersen<dimension, float_type>::verlet_nvt_andersen(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type const> box
  , shared_ptr<random_type> random
  , float_type timestep
  , float_type temperature
  , float_type coll_rate
  , shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , random_(random)
  , coll_rate_(coll_rate)
  , logger_(logger)
{
    this->timestep(timestep);
    this->temperature(temperature);
    LOG("collision rate with heat bath: " << coll_rate_);
}

template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::timestep(double timestep)
{
    timestep_ = static_cast<float_type>(timestep);
    timestep_half_ = 0.5 * timestep;
    coll_prob_ = coll_rate_ * timestep;

    LOG("integration timestep: " << timestep_);
}

template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::temperature(double temperature)
{
    temperature_ = static_cast<float_type>(temperature);
    sqrt_temperature_ = sqrt(temperature_);

    LOG("temperature of heat bath: " << temperature_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::integrate()
{
    scoped_timer_type timer(runtime_.integrate);

    for (size_t i = 0; i < particle_->nbox; ++i) {
        unsigned int type = particle_->type[i];
        float_type mass = particle_->mass[type];
        vector_type& v = particle_->v[i] += particle_->f[i] * timestep_half_ / mass;
        vector_type& r = particle_->r[i] += v * timestep_;
        // enforce periodic boundary conditions
        // TODO: reduction is now to (-L/2, L/2) instead of (0, L) as before
        // check that this is OK
        particle_->image[i] += box_->reduce_periodic(r);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::finalize()
{
    scoped_timer_type timer(runtime_.finalize);

    // cache random numbers
    float_type rng_cache = 0;
    bool rng_cache_valid = false;

    // loop over all particles
    for (size_t i = 0; i < particle_->nbox; ++i) {
        // is deterministic step?
        if (random_->uniform<float_type>() > coll_prob_) {
            unsigned int type = particle_->type[i];
            float_type mass = particle_->mass[type];
            particle_->v[i] += particle_->f[i] * timestep_half_ / mass;
        }
        // stochastic coupling with heat bath
        else {
            // assign two velocity components at a time
            vector_type& v = particle_->v[i];
            for (unsigned i=0; i < dimension-1; i+=2) {
                tie(v[i], v[i+1]) = random_->normal(sqrt_temperature_);
            }
            // handle last component separately for odd dimensions
            if (dimension % 2 == 1) {
                if (rng_cache_valid) {
                    v[dimension-1] = rng_cache;
                }
                else {
                    tie(v[dimension-1], rng_cache) = random_->normal(sqrt_temperature_);
                }
                rng_cache_valid = !rng_cache_valid;
            }
        }
    }
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(verlet_nvt_andersen<dimension, float_type> const&)
{
    return verlet_nvt_andersen<dimension, float_type>::module_name();
}

template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::luaopen(lua_State* L)
{
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
                    class_<verlet_nvt_andersen, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type>
                          , shared_ptr<box_type const>
                          , shared_ptr<random_type>
                          , float_type
                          , float_type
                          , float_type
                          , shared_ptr<logger_type>
                        >())
                        .property("collision_rate", &verlet_nvt_andersen::collision_rate)
                        .property("module_name", &module_name_wrapper<dimension, float_type>)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("integrate", &runtime::integrate)
                                .def_readonly("finalize", &runtime::finalize)
                        ]
                        .def_readonly("runtime", &verlet_nvt_andersen::runtime_)
                ]
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

} // namespace mdsim
} // namespace host
} // namespace integrators
} // namespace halmd
