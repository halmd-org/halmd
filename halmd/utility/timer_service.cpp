/*
 * Copyright © 2011  Peter Colberg
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

#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/timer_service.hpp>

namespace halmd {
namespace utility {

connection timer_service::on_periodic(slot_function_type const& slot, time_type interval)
{
    return on_periodic(slot, interval, interval);
}

connection timer_service::on_periodic(slot_function_type const& slot, time_type interval, time_type start)
{
    event e;
    e.time_ = time(0) + start;
    e.interval_ = interval;
    e.slot_ = slot;
    return on_periodic_.connect(e);
}

void timer_service::process()
{
    on_periodic_(time(0));
}

void timer_service::event::operator()(time_type const& time)
{
    if (time >= time_) {
        slot_();
        time_ += interval_;
    }
}

void timer_service::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("utility")
        [
            class_<timer_service, std::shared_ptr<timer_service> >("timer_service")
                .def(constructor<>())
                .def("on_periodic", static_cast<connection (timer_service::*)(slot_function_type const&, time_type)>(&timer_service::on_periodic))
                .def("on_periodic", static_cast<connection (timer_service::*)(slot_function_type const&, time_type, time_type)>(&timer_service::on_periodic))
                .def("process", &timer_service::process)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_utility_timer_service(lua_State* L)
{
    timer_service::luaopen(L);
    return 0;
}

} // namespace utility
} // namespace halmd
