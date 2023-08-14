/*
 * Copyright © 2017 Arthur Straube
 * Copyright © 2017 Felix Höfling
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

#include <halmd/mdsim/gpu/positions/random.hpp>
#include <halmd/mdsim/gpu/positions/random_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace positions {

template <int dimension, typename float_type, typename RandomNumberGenerator>
random<dimension, float_type, RandomNumberGenerator>::random(
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

template <int dimension, typename float_type, typename RandomNumberGenerator>
void random<dimension, float_type, RandomNumberGenerator>::set()
{
    auto position = make_cache_mutable(particle_->position());
    auto image = make_cache_mutable(particle_->image());

    LOG_TRACE("randomly distributing positions of " << position->size() << " particles");

    scoped_timer_type timer(runtime_.set);

    // edge lengths of cuboid slab centred around the origin
    vector_type slab_length = element_prod(static_cast<vector_type>(box_->length()), slab_);

    try {
        auto& random_kernel = random_wrapper<dimension, typename rng_type::rng_type>::kernel.uniform;

        random_kernel.configure(rng_->rng().dim.grid, rng_->rng().dim.block);
        random_kernel(
            &*position->begin()
          , particle_->nparticle()
          , particle_->dim().threads()
          , slab_length
          , rng_->rng().rng()
        );

        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to generate random particle positions on GPU");
        throw;
    }

    // reset particle image vectors
    cuda::memset(image->begin(), image->begin() + image->capacity(), 0);
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void random<dimension, float_type, RandomNumberGenerator>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("positions")
            [
                class_<random, _Base>()
                    .property("slab", &random::slab)
                    .def("set", &random::set)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("set", &runtime::set)
                    ]
                    .def_readonly("runtime", &random::runtime_)
              , def("random", &std::make_shared<random
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

using halmd::random::gpu::rand48;

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_positions_random(lua_State* L)
{
    random<3, float, rand48>::luaopen(L);
    random<2, float, rand48>::luaopen(L);
    return 0;
}

// explicit instantiation
template class random<3, float, rand48>;
template class random<2, float, rand48>;

} // namespace positions
} // namespace gpu
} // namespace mdsim
} // namespace halmd
