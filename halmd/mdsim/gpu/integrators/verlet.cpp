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
#include <halmd/mdsim/gpu/integrators/verlet.hpp>
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
verlet<dimension, float_type>::verlet(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , double timestep
)
  // dependency injection
  : particle(particle)
  , box(box)
  // reference CUDA C++ verlet_wrapper
  , wrapper(&verlet_wrapper<dimension>::wrapper)
{
    this->timestep(timestep);

    try {
        cuda::copy(static_cast<vector_type>(box->length()), wrapper->box_length);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to initialize Verlet integrator symbols");
        throw;
    }
}

/**
 * set integration time-step
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::timestep(double timestep)
{
    timestep_ = timestep;
    timestep_half_ = 0.5 * timestep_;

    try {
        cuda::copy(timestep_, wrapper->timestep);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to initialize Verlet integrator symbols");
        throw;
    }

    LOG("integration timestep: " << timestep_);
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.integrate, "integrate", "first half-step of velocity-Verlet");
    profiler.register_runtime(runtime_.finalize, "finalize", "second half-step of velocity-Verlet");
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::integrate()
{
    try {
        scoped_timer<timer> timer_(runtime_.integrate);
        cuda::configure(particle->dim.grid, particle->dim.block);
        wrapper->integrate(
            particle->g_r, particle->g_image, particle->g_v, particle->g_f);
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
        scoped_timer<timer> timer_(runtime_.finalize);
        cuda::configure(particle->dim.grid, particle->dim.block);
        wrapper->finalize(particle->g_v, particle->g_f);
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream second leapfrog step on GPU");
        throw;
    }
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(verlet<dimension, float_type> const&)
{
    return verlet<dimension, float_type>::module_name();
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("integrators")
                [
                    class_<verlet, shared_ptr<_Base>, bases<_Base> >(class_name.c_str())
                        .def(constructor<shared_ptr<particle_type>, shared_ptr<box_type>, double>())
                        .def("register_runtimes", &verlet::register_runtimes)
                        .property("module_name", &module_name_wrapper<dimension, float_type>)
                ]
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
