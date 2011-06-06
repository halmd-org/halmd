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

#include <cmath>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Initialize simulation
 */
template <int dimension>
core<dimension>::core()
  // dependency injection
  : clock(make_shared<clock_type>())
{
    LOG("dimension of positional coordinates: " << dimension);
}

/**
 * register module runtime accumulators
 */
template <int dimension>
void core<dimension>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.prepare, "prepare", "microscopic state preparation");
    profiler.register_runtime(runtime_.mdstep, "mdstep", "MD integration step");
}

/**
 * Prepare microscopic system state
 */
template <int dimension>
void core<dimension>::prepare()
{
    scoped_timer<timer> timer_(runtime_.prepare);
    particle->set();
    position->set();
    velocity->set();
    if (neighbour) {
        if (sort) {
            sort->order();
        }
        neighbour->update();
    }
    force->compute();
}

/**
 * Perform a single MD integration step
 */
template <int dimension>
void core<dimension>::mdstep()
{
    scoped_timer<timer> timer_(runtime_.mdstep);

    LOG_TRACE("performing MD step #" << clock->step() + 1); //< output 1-based counter consistent with output files

    integrator->integrate();
    if (neighbour && neighbour->check()) {
        if (sort) {
            sort->order();
        }
        neighbour->update();
    }
    force->compute();
    integrator->finalize();

    clock->advance(integrator->timestep()); // FIXME move to beginning of function?
}

/**
 * wrap dimension template parameter for Lua
 */
template <int dimension>
static int get_dimension(core<dimension> const&)
{
    return dimension;
}

template <int dimension>
void core<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("core_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<core, shared_ptr<core> >(class_name.c_str())
                .def(constructor<>())
                .def("register_runtimes", &core::register_runtimes)
                .def_readwrite("particle", &core::particle)
                .def_readwrite("box", &core::box)
                .def_readwrite("force", &core::force)
                .def_readwrite("neighbour", &core::neighbour)
                .def_readwrite("sort", &core::sort)
                .def_readwrite("integrator", &core::integrator)
                .def_readwrite("position", &core::position)
                .def_readwrite("velocity", &core::velocity)
                .def_readonly("clock", &core::clock)
                .property("dimension", &get_dimension<dimension>)
                .property("step", &core::step)
                .def("prepare", &core::prepare)
                .def("mdstep", &core::mdstep)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_core(lua_State* L)
{
    core<3>::luaopen(L);
    core<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class core<3>;
template class core<2>;

} // namespace mdsim

} // namespace halmd