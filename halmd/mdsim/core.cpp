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

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {

/**
 * Initialize simulation
 */
core::core(boost::shared_ptr<clock_type> clock)
  // dependency injection
  : clock_(clock)
{
}

/**
 * Prepare microscopic system state
 */
void core::setup()
{
    scoped_timer_type timer(runtime_.setup);

    on_prepend_setup_();
    on_setup_();
    on_append_setup_();
}

/**
 * Perform a single MD integration step
 */
void core::mdstep()
{
    scoped_timer_type timer(runtime_.mdstep);

    // increment 1-based simulation step
    clock_->advance();

    LOG_TRACE("performing MD step #" << clock_->step());

    on_prepend_integrate_();
    on_integrate_();
    on_append_integrate_();
    on_prepend_force_();
    on_force_();
    on_append_force_();
    on_prepend_finalize_();
    on_finalize_();
    on_append_finalize_();
}

void core::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<core, boost::shared_ptr<core> >("core")
                .def(constructor<boost::shared_ptr<core::clock_type> >())
                .def("setup", &core::setup)
                .def("mdstep", &core::mdstep)
                .def("on_prepend_setup", &core::on_prepend_setup)
                .def("on_setup", &core::on_setup)
                .def("on_append_setup", &core::on_append_setup)
                .def("on_prepend_integrate", &core::on_prepend_integrate)
                .def("on_integrate", &core::on_integrate)
                .def("on_append_integrate", &core::on_append_integrate)
                .def("on_prepend_force", &core::on_prepend_force)
                .def("on_force", &core::on_force)
                .def("on_append_force", &core::on_append_force)
                .def("on_prepend_finalize", &core::on_prepend_finalize)
                .def("on_finalize", &core::on_finalize)
                .def("on_append_finalize", &core::on_append_finalize)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("setup", &runtime::setup)
                        .def_readonly("mdstep", &runtime::mdstep)
                ]
                .def_readonly("runtime", &core::runtime_)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_core(lua_State* L)
{
    core::luaopen(L);
    return 0;
}

} // namespace mdsim
} // namespace halmd
