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

#include <halmd/mdsim/core.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {

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

static std::function<void ()>
wrap_setup(boost::shared_ptr<core> self)
{
    return [=]() {
        self->setup();
    };
}

static std::function<void ()>
wrap_mdstep(boost::shared_ptr<core> self)
{
    return [=]() {
        self->mdstep();
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_prepend_setup(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_prepend_setup(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_setup(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_setup(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_append_setup(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_append_setup(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_prepend_integrate(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_prepend_integrate(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_integrate(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_integrate(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_append_integrate(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_append_integrate(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_prepend_force(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_prepend_force(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_force(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_force(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_append_force(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_append_force(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_prepend_finalize(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_prepend_finalize(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_finalize(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_finalize(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_append_finalize(boost::shared_ptr<core> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_append_finalize(slot);
    };
}

void core::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<core, boost::shared_ptr<core> >("core")
                .def(constructor<>())
                .property("setup", &wrap_setup)
                .property("mdstep", &wrap_mdstep)
                .property("on_prepend_setup", &wrap_on_prepend_setup)
                .property("on_setup", &wrap_on_setup)
                .property("on_append_setup", &wrap_on_append_setup)
                .property("on_prepend_integrate", &wrap_on_prepend_integrate)
                .property("on_integrate", &wrap_on_integrate)
                .property("on_append_integrate", &wrap_on_append_integrate)
                .property("on_prepend_force", &wrap_on_prepend_force)
                .property("on_force", &wrap_on_force)
                .property("on_append_force", &wrap_on_append_force)
                .property("on_prepend_finalize", &wrap_on_prepend_finalize)
                .property("on_finalize", &wrap_on_finalize)
                .property("on_append_finalize", &wrap_on_append_finalize)
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
