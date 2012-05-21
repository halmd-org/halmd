/*
 * Copyright © 2011-2012  Michael Kopp and Felix Höfling
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

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/integrators/euler.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type>
euler<dimension, float_type>::euler(
    boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<box_type const> box
  , double timestep
  , boost::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
{
    set_timestep(timestep);
}

/**
 * set integration timestep
 */
template <int dimension, typename float_type>
void euler<dimension, float_type>::set_timestep(double timestep)
{
    timestep_ = timestep;
}

/**
 * perform Euler integration: update positions from velocities
 */
template <int dimension, typename float_type>
void euler<dimension, float_type>::integrate()
{
    try {
        scoped_timer_type timer(runtime_.integrate);
        cuda::configure(particle_->dim.grid, particle_->dim.block);
        wrapper_type::kernel.integrate(
            particle_->position(), particle_->image(), particle_->velocity()
          , timestep_
          , static_cast<vector_type>(box_->length())
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream Euler integration on GPU");
        throw;
    }
}

template <typename integrator_type>
static boost::function<void ()>
wrap_integrate(boost::shared_ptr<integrator_type> self)
{
    return bind(&integrator_type::integrate, self);
}

template <typename integrator_type>
static boost::function<void ()>
wrap_finalize(boost::shared_ptr<integrator_type> self)
{
    return bind(&integrator_type::finalize, self);
}

template <int dimension, typename float_type>
void euler<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<euler>()
                    .property("integrate", &wrap_integrate<euler>)
                    .property("timestep", &euler::timestep)
                    .def("set_timestep", &euler::set_timestep)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("integrate", &runtime::integrate)
                    ]
                    .def_readonly("runtime", &euler::runtime_)

              , def("euler", &boost::make_shared<euler
                  , boost::shared_ptr<particle_type>
                  , boost::shared_ptr<box_type const>
                  , double
                  , boost::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_integrators_euler(lua_State* L)
{
    euler<3, float>::luaopen(L);
    euler<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class euler<3, float>;
template class euler<2, float>;

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
