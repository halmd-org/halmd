/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/integrators/verlet_nvt_andersen.hpp>
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
verlet_nvt_andersen<dimension, float_type>::verlet_nvt_andersen(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<random_type> random
  , float_type timestep, float_type temperature, float_type coll_rate
)
  // dependency injection
  : particle(particle)
  , box(box)
  , random(random)
  , coll_rate_(coll_rate)
{
    this->timestep(timestep);
    this->temperature(temperature);
    LOG("collision rate with heat bath: " << coll_rate_);
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void verlet_nvt_andersen<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
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
    scoped_timer<timer> timer_(at_key<integrate_>(runtime_));

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
void verlet_nvt_andersen<dimension, float_type>::finalize()
{
    scoped_timer<timer> timer_(at_key<finalize_>(runtime_));

    // cache random numbers
    float_type rng_cache = 0;
    bool rng_cache_valid = false;

    // loop over all particles
    for (size_t i = 0; i < particle->nbox; ++i) {
        // is deterministic step?
        if (random->uniform<float_type>() > coll_prob_) {
            particle->v[i] += particle->f[i] * timestep_half_;
        }
        // stochastic coupling with heat bath
        else {
            // assign two velocity components at a time
            vector_type& v = particle->v[i];
            for (unsigned i=0; i < dimension-1; i+=2) {
                tie(v[i], v[i+1]) = random->normal(sqrt_temperature_);
            }
            // handle last component separately for odd dimensions
            if (dimension % 2 == 1) {
                if (rng_cache_valid) {
                    v[dimension-1] = rng_cache;
                }
                else {
                    tie(v[dimension-1], rng_cache) = random->normal(sqrt_temperature_);
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
                            verlet_nvt_andersen
                          , shared_ptr<_Base_Base>
                          , bases<_Base_Base, _Base>
                        >(class_name.c_str())
                            .def(constructor<
                                shared_ptr<particle_type>
                              , shared_ptr<box_type>
                              , shared_ptr<random_type>
                              , float_type, float_type, float_type>()
                            )
                            .def("register_runtimes", &verlet_nvt_andersen::register_runtimes)
                            .property("collision_rate", &verlet_nvt_andersen::collision_rate)
                            .property("module_name", &module_name_wrapper<dimension, float_type>)
                    ]
                ]
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &verlet_nvt_andersen<3, double>::luaopen
    ]
    [
        &verlet_nvt_andersen<2, double>::luaopen
    ];
#else
    [
        &verlet_nvt_andersen<3, float>::luaopen
    ]
    [
        &verlet_nvt_andersen<2, float>::luaopen
    ];
#endif
}

} // namespace

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class verlet_nvt_andersen<3, double>;
template class verlet_nvt_andersen<2, double>;
#else
template class verlet_nvt_andersen<3, float>;
template class verlet_nvt_andersen<2, float>;
#endif

}}} // namespace mdsim::host::integrators

} // namespace halmd
