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
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <cmath>

#include <halmd/mdsim/gpu/integrators/verlet.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type>
verlet<dimension, float_type>::verlet(
    boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<box_type const> box
  , double timestep
  , boost::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
  // reference CUDA C++ verlet_wrapper
  , wrapper_(&verlet_wrapper<dimension>::wrapper)
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
    try {
        scoped_timer_type timer(runtime_.integrate);
        cuda::configure(
            particle_->dim.grid, particle_->dim.block
        );
        wrapper_->integrate(
            particle_->position(), particle_->image(), particle_->velocity()
          , particle_->force()
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
    // TODO: possibly a performance critical issue:
    // the old implementation had this loop included in update_forces(),
    // which saves one additional read of the forces plus the additional kernel execution
    // and scheduling
    try {
        scoped_timer_type timer(runtime_.finalize);
        cuda::configure(
            particle_->dim.grid, particle_->dim.block
        );
        wrapper_->finalize(
            particle_->velocity()
          , particle_->force()
          , timestep_
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream second leapfrog step on GPU");
        throw;
    }
}

template <typename integrator_type>
static std::function<void ()>
wrap_integrate(boost::shared_ptr<integrator_type> self)
{
    return bind(&integrator_type::integrate, self);
}

template <typename integrator_type>
static std::function<void ()>
wrap_finalize(boost::shared_ptr<integrator_type> self)
{
    return bind(&integrator_type::finalize, self);
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<verlet>()
                    .property("integrate", &wrap_integrate<verlet>)
                    .property("finalize", &wrap_finalize<verlet>)
                    .property("timestep", &verlet::timestep)
                    .def("set_timestep", &verlet::set_timestep)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("integrate", &runtime::integrate)
                            .def_readonly("finalize", &runtime::finalize)
                    ]
                    .def_readonly("runtime", &verlet::runtime_)

              , def("verlet", &boost::make_shared<verlet
                  , boost::shared_ptr<particle_type>
                  , boost::shared_ptr<box_type const>
                  , double
                  , boost::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_integrators_verlet(lua_State* L)
{
    verlet<3, float>::luaopen(L);
    verlet<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class verlet<3, float>;
template class verlet<2, float>;

} // namespace mdsim
} // namespace gpu
} // namespace integrators
} // namespace halmd
