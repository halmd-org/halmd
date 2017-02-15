/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/gpu/integrators/verlet.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type>
verlet<dimension, float_type>::verlet(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , double timestep
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
  // reference CUDA C++ verlet_wrapper
  , wrapper_(&verlet_wrapper<dimension, float_type>::wrapper)
{
    set_timestep(timestep);
}

/**
 * set integration time-step
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::set_timestep(double timestep)
{
    timestep_ = timestep;
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::integrate()
{
    LOG_TRACE("update positions and velocities");

    force_array_type const& force = read_cache(particle_->force());

    // invalidate the particle caches after accessing the force!
    auto position = make_cache_mutable(particle_->position());
    auto velocity = make_cache_mutable(particle_->velocity());
    auto image = make_cache_mutable(particle_->image());

    scoped_timer_type timer(runtime_.integrate);

    try {
        cuda::configure(particle_->dim().grid, particle_->dim().block);
        wrapper_->integrate(
            position->data()
          , image->data()
          , velocity->data()
          , force.data()
          , timestep_
          , static_cast<vector_type>(box_->length())
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream first leapfrog step on GPU");
        throw;
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::finalize()
{
    LOG_TRACE("update velocities");

    force_array_type const& force = read_cache(particle_->force());

    // invalidate the particle caches after accessing the force!
    auto velocity = make_cache_mutable(particle_->velocity());

    scoped_timer_type timer(runtime_.finalize);

    try {
        cuda::configure(particle_->dim().grid, particle_->dim().block);
        wrapper_->finalize(
            velocity->data()
          , force.data()
          , timestep_
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream second leapfrog step on GPU");
        throw;
    }
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<verlet>()
                    .def("integrate", &verlet::integrate)
                    .def("finalize", &verlet::finalize)
                    .def("set_timestep", &verlet::set_timestep)
                    .property("timestep", &verlet::timestep)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("integrate", &runtime::integrate)
                            .def_readonly("finalize", &runtime::finalize)
                    ]
                    .def_readonly("runtime", &verlet::runtime_)

              , def("verlet", &std::make_shared<verlet
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type const>
                  , double
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_integrators_verlet(lua_State* L)
{
    verlet<3, float>::luaopen(L);
    verlet<2, float>::luaopen(L);
    verlet<3, dsfloat>::luaopen(L);
    verlet<2, dsfloat>::luaopen(L);
    return 0;
}

// explicit instantiation
template class verlet<3, float>;
template class verlet<2, float>;
template class verlet<3, dsfloat>;
template class verlet<2, dsfloat>;

} // namespace mdsim
} // namespace gpu
} // namespace integrators
} // namespace halmd
