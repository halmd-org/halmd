/*
 * Copyright © 2011-2022 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#include <halmd/observables/host/kinetic_energy_density_mode.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
kinetic_energy_density_mode<dimension, float_type>::kinetic_energy_density_mode(
    shared_ptr<particle_type const> particle
  , shared_ptr<particle_group_type> particle_group
  , shared_ptr<wavevector_type const> wavevector
  , shared_ptr<logger> logger
)
    // dependency injection
  : particle_(particle)
  , particle_group_(particle_group)
  , wavevector_(wavevector)
  , logger_(logger)
{}

/**
 * Acquire density modes from particle positions
 */
template <int dimension, typename float_type>
shared_ptr<typename kinetic_energy_density_mode<dimension, float_type>::result_type const>
kinetic_energy_density_mode<dimension, float_type>::acquire()
{
    // check validity of caches
    auto const& group_cache  = particle_group_->ordered();
    auto const& position_cache = particle_->position();
    auto const& velocity_cache = particle_->velocity();
    auto const& mass_cache = particle_->mass();

    if (group_cache_ != group_cache || position_cache_ != position_cache
     || velocity_cache_ != velocity_cache || mass_cache_ != mass_cache
    ) {
        // obtain read access to input caches
        auto const& group = read_cache(group_cache);
        auto const& position = read_cache(position_cache);
        auto const& velocity = read_cache(velocity_cache);
        auto const& mass = read_cache(mass_cache);

        LOG_TRACE("acquire sample");

        scoped_timer_type timer(runtime_.acquire);

        auto const& wavevector = wavevector_->value(); // array of wavevectors

        // allocate new memory which allows modules (e.g.,
        // dynamics::blocking_scheme) to hold a previous copy of the result or
        // to track the update via std::weak_ptr.
        result_ = make_shared<result_type>(wavevector.size());

        // compute density modes
        // initialise result array
        fill(begin(*result_), end(*result_), 0);

        // compute sum of exponentials: rho_q = sum_r exp(-i q·r)
        // 1st loop: iterate over particle group
        for (auto i : group) {
            vector_type const& r = position[i];
            vector_type const& v = velocity[i];
            float_type m = mass[i];
            // 2nd loop: iterate over wavevectors
            typedef typename result_type::value_type mode_type;
            auto rho_q = begin(*result_);
            for (auto const& q : wavevector) {
                float_type q_r = inner_prod(static_cast<vector_type>(q), r);
                float_type ekin = m * inner_prod(v, v) / 2;
                *rho_q++ += mode_type({{ ekin * cos(q_r), -ekin * sin(q_r) }});
            }
        }

        // update cache observers
        group_cache_ = group_cache;
        position_cache_ = position_cache;
        velocity_cache_ = velocity_cache;
        mass_cache_ = mass_cache;
   }

    return result_;
}

template <int dimension, typename float_type>
void kinetic_energy_density_mode<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<kinetic_energy_density_mode>()
                    .property("acquisitor", &kinetic_energy_density_mode::acquisitor)
                    .property("wavevector", &kinetic_energy_density_mode::wavevector)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("acquire", &runtime::acquire)
                    ]
                    .def_readonly("runtime", &kinetic_energy_density_mode::runtime_)
            ]
          , def("kinetic_energy_density_mode", &make_shared<kinetic_energy_density_mode
              , shared_ptr<particle_type const>
              , shared_ptr<particle_group_type>
              , shared_ptr<wavevector_type const>
              , shared_ptr<logger>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_kinetic_energy_density_mode(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    kinetic_energy_density_mode<3, double>::luaopen(L);
    kinetic_energy_density_mode<2, double>::luaopen(L);
#else
    kinetic_energy_density_mode<3, float>::luaopen(L);
    kinetic_energy_density_mode<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class kinetic_energy_density_mode<3, double>;
template class kinetic_energy_density_mode<2, double>;
#else
template class kinetic_energy_density_mode<3, float>;
template class kinetic_energy_density_mode<2, float>;
#endif

}  // namespace host
}  // namespace observables
}  // namespace halmd
