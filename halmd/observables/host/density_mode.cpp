/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#include <halmd/observables/host/density_mode.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
density_mode<dimension, float_type>::density_mode(
    shared_ptr<particle_type const> particle
  , shared_ptr<particle_group_type> particle_group
  , shared_ptr<wavevector_type const> wavevector
  , shared_ptr<logger_type> logger
)
    // dependency injection
  : particle_(particle)
  , particle_group_(particle_group)
  , wavevector_(wavevector)
  , logger_(logger)
{
}

/**
 * Acquire sample of all density modes from particle group
 */
template <int dimension, typename float_type>
shared_ptr<typename density_mode<dimension, float_type>::sample_type const>
density_mode<dimension, float_type>::acquire()
{
    scoped_timer_type timer(runtime_.acquire);

    LOG_TRACE("acquire sample");

    auto const& wavevector = wavevector_->value(); // array of wavevectors

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    rho_sample_ = make_shared<sample_type>(wavevector.size());

    // compute density modes
    mode_array_type& rho_vector = rho_sample_->rho();
    // initialise result array
    fill(begin(rho_vector), end(rho_vector), 0);

    // compute sum of exponentials: rho_q = sum_r exp(-i q·r)
    // 1st loop: iterate over particles
    auto const& unordered = read_cache(particle_group_->unordered());
    auto const& position  = read_cache(particle_->position());

    for (typename particle_group_type::size_type i : unordered) {
        vector_type const& r = position[i];
        // 2nd loop: iterate over wavevectors
        auto rho_q = begin(rho_vector);
        for (auto const& q : wavevector) {
            float_type q_r = inner_prod(q, r);
            *rho_q++ += mode_type({{ cos(q_r), -sin(q_r) }});
        }
    }

    return rho_sample_;
}

template <typename sample_type, typename density_mode_type>
static function<shared_ptr<sample_type const> ()>
wrap_acquire(shared_ptr<density_mode_type> density_mode)
{
    return [=]() {
       return density_mode->acquire();
    };
}

template <int dimension, typename float_type>
void density_mode<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<density_mode>()
                    .def("acquire", &wrap_acquire<sample_type, density_mode>)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("acquire", &runtime::acquire)
                    ]
                    .def_readonly("runtime", &density_mode::runtime_)
            ]
          , def("density_mode", &make_shared<density_mode
              , shared_ptr<particle_type const>
              , shared_ptr<particle_group_type>
              , shared_ptr<wavevector_type const>
              , shared_ptr<logger_type>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_density_mode(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    density_mode<3, double>::luaopen(L);
    density_mode<2, double>::luaopen(L);
#else
    density_mode<3, float>::luaopen(L);
    density_mode<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class density_mode<3, double>;
template class density_mode<2, double>;
#else
template class density_mode<3, float>;
template class density_mode<2, float>;
#endif

}  // namespace host
}  // namespace observables
}  // namespace halmd
