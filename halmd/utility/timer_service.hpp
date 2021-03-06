/*
 * Copyright © 2011  Peter Colberg
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

#ifndef HALMD_UTILITY_TIMER_SERVICE_HPP
#define HALMD_UTILITY_TIMER_SERVICE_HPP

#include <ctime>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace utility {

/**
 * This class provides a periodic timer service.
 *
 * The timer service may be used during the simulation to perform periodic
 * tasks such as flushing HDF5 buffers to disk, writing a snapshot of the
 * simulation state to disk, or updating the runtime estimator.
 */
class timer_service
{
public:
    typedef std::time_t time_type;
    typedef signal<void ()>::slot_function_type slot_function_type;

    /** connect slot to be invoked periodically */
    connection on_periodic(slot_function_type const& slot, time_type interval);
    /** connect slot to be invoked periodically, starting at given relative time */
    connection on_periodic(slot_function_type const& slot, time_type interval, time_type start);
    /** process timer event queue */
    void process();
    /** Lua bindings */
    static void luaopen(lua_State* L);

private:
    /**
     * Periodic timed event
     */
    struct event
    {
        /** process event */
        void operator()(time_type const& time);
        /** time in seconds since the epoch */
        time_type time_;
        /** time interval in seconds */
        time_type interval_;
        /** slot */
        slot_function_type slot_;
    };

    signal<void (time_type const&)> on_periodic_;
};

} // namespace utility
} // namespace halmd

#endif /* ! HALMD_UTILITY_TIMER_SERVICE_HPP */
