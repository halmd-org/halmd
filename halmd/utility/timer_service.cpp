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

/**
 * Register callback to be invoked periodically.
 *
 * @param interval periodic interval in seconds
 * @param slot callback function or functor
 */
void timer_service::on_periodic(slot_function_type const& slot, time_t interval)
{
    event e;
    e.time = time(0) + interval;
    e.interval = interval;
    e.slot = slot;
    queue_.push(e);
}

/**
 * Process timer event queue.
 */
void timer_service::process()
{
    time_t const t = time(0);
    while (!queue_.empty() && t >= queue_.top().time) {
        event e(queue_.top());
        e.slot();
        queue_.pop();
        e.time += e.interval;
        queue_.push(e);
    }
}

signal<void (uint64_t)>::slot_function_type
process_wrapper(shared_ptr<timer_service> ts)
{
    return bind(&timer_service::process, ts);
}

HALMD_LUA_API int luaopen_libhalmd_utility_timer_service(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("utility")
        [
            class_<timer_service, shared_ptr<timer_service> >("timer_service")
                .def(constructor<>())
                .def("on_periodic", &timer_service::on_periodic)
                .property("process", &process_wrapper)
        ]
    ];
    return 0;
}

} // namespace utility
} // namespace halmd
