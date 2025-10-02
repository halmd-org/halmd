/*
 * Copyright © 2016 Arthur Straube
 * Copyright © 2016 Felix Höfling
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
#include <functional>

#include <halmd/mdsim/host/positions/uniform.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace positions {

template <int dimension, typename float_type>
uniform<dimension, float_type>::uniform(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<rng_type> rng
  , vector_type const& slab
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , rng_(rng)
  , logger_(logger)
  , slab_(slab)
{
    // FIXME replace slab by future 'geometry' modules from Nicolas
    if (*min_element(slab_.begin(), slab_.end()) <= 0 ||
        *max_element(slab_.begin(), slab_.end()) > 1
       ) {
        throw std::logic_error("slab extents must be a fraction between 0 and 1");
    }

    if (*min_element(slab_.begin(), slab_.end()) < 1) {
        LOG("restrict initial particle positions to slab: " << slab_);
    }
}

template <int dimension, typename float_type>
void uniform<dimension, float_type>::set()
{
    auto position = make_cache_mutable(particle_->position());
    auto image = make_cache_mutable(particle_->image());

    LOG_TRACE("uniformly distributing random positions of " << position->size() << " particles");

    scoped_timer_type timer(runtime_.set);

    // edge lengths of cuboid slab centred around the origin
    vector_type length = element_prod(box_->length(), slab_);

    // iterate over all particles
    for (auto &r : *position) {
        // assign to each component uniform random values from [-1/2, 1/2)
        for (unsigned int i = 0; i < dimension; ++i) {
            r[i] = rng_->uniform<float_type>() - float_type(.5);
        }
        // scale each component by slab size
        r = element_prod(r, length);
    }

    // reset particle image vectors
    fill(image->begin(), image->end(), 0);
}

template <int dimension, typename float_type>
void uniform<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("positions")
            [
                class_<uniform, _Base>()
                    .property("slab", &uniform::slab)
                    .def("set", &uniform::set)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("set", &runtime::set)
                    ]
                    .def_readonly("runtime", &uniform::runtime_)
              , def("uniform", &std::make_shared<uniform
                  , std::shared_ptr<particle_type>
                    , std::shared_ptr<box_type const>
                    , std::shared_ptr<rng_type>
                    , vector_type const&
                    , std::shared_ptr<logger>
                  >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_positions_uniform(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    uniform<3, double>::luaopen(L);
    uniform<2, double>::luaopen(L);
#else
    uniform<3, float>::luaopen(L);
    uniform<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class uniform<3, double>;
template class uniform<2, double>;
#else
template class uniform<3, float>;
template class uniform<2, float>;
#endif

} // namespace positions
} // namespace host
} // namespace mdsim
} // namespace halmd
