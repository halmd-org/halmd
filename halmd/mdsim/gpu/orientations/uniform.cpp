/*
 * Copyright © 2017-2018  Jake Atwell
 * Copyright © 2016       Manuel Dibak
 * Copyright © 2008-2011  Felix Höfling
 * Copyright © 2013       Nicolas Höft
 * Copyright © 2008-2011  Peter Colberg
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

#include <halmd/mdsim/gpu/orientations/uniform_kernel.hpp>
#include <halmd/mdsim/gpu/orientations/uniform.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace orientations {

template <int dimension, typename float_type, typename RandomNumberGenerator>
uniform<dimension, float_type, RandomNumberGenerator>::uniform(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<random_type> random
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , random_(random)
  , logger_(logger)
{}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void uniform<dimension, float_type, RandomNumberGenerator>::set()
{
    scoped_timer_type timer(runtime_.set);

    auto orientation = make_cache_mutable(particle_->orientation());
    size_t npart = particle_->nparticle();

    try {
        cuda::configure(random_->rng().dim.grid, random_->rng().dim.block);
        uniform_wrapper<dimension, float_type, rng_type>::kernel.uniform(
            orientation->data(), npart, random_->rng().rng()
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to generate uniform particle orientation on GPU");
        throw;
    }
    LOG("setting particle orientations");
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void uniform<dimension, float_type, RandomNumberGenerator>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("orientations")
            [
                class_<uniform>()
                    .def("set", &uniform::set)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("set", &runtime::set)
                    ]
                    .def_readonly("runtime", &uniform::runtime_)
              , def("uniform", &std::make_shared<uniform
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<random_type>
                  >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_orientations_uniform(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    uniform<3, float, halmd::random::gpu::rand48>::luaopen(L);
    uniform<2, float, halmd::random::gpu::rand48>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    uniform<3, dsfloat, halmd::random::gpu::rand48>::luaopen(L);
    uniform<2, dsfloat, halmd::random::gpu::rand48>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class uniform<3, float, halmd::random::gpu::rand48>;
template class uniform<2, float, halmd::random::gpu::rand48>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class uniform<3, dsfloat, halmd::random::gpu::rand48>;
template class uniform<2, dsfloat, halmd::random::gpu::rand48>;
#endif

} // namespace mdsim
} // namespace gpu
} // namespace orientations
} // namespace halmd
