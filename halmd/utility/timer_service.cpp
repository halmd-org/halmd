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

#include <boost/bind.hpp>

#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/timer_service.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace utility {

timer_service::connection_type
timer_service::on_periodic(slot_function_type const& slot, time_type interval)
{
    event e;
    e.time_ = time(0) + interval;
    e.interval_ = interval;
    e.slot_ = slot;
    return on_periodic_.connect(e);
}

void timer_service::process()
{
    on_periodic_();
}

void timer_service::event::operator()()
{
    if (time(0) >= time_) {
        slot_();
        time_ += interval_;
    }
}

signal<void (uint64_t)>::slot_function_type
wrap_process(shared_ptr<timer_service> ts)
{
    return bind(&timer_service::process, ts);
}

void timer_service::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("utility")
        [
            class_<timer_service, shared_ptr<timer_service> >("timer_service")
                .def(constructor<>())
                .def("on_periodic", &timer_service::on_periodic)
                .property("process", &wrap_process)
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
