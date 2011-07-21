/*
 * Copyright © 2010-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_UTILITY_PROFILER_HPP
#define HALMD_UTILITY_PROFILER_HPP

#include <lua.hpp>
#include <vector>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd {
namespace utility {

class profiler
{
private:
    typedef signal<void ()> signal_type;

public:
    typedef timer timer_type;
    typedef accumulator<double> accumulator_type;
    typedef scoped_timer<timer_type> scoped_timer_type;
    typedef signal_type::slot_function_type slot_function_type;
    typedef signal_type::connection connection_type;

    /** logs timer resolution */
    profiler();
    /** connect accumulator for profiling */
    connection_type on_profile(boost::shared_ptr<accumulator_type> acc, std::string const& desc);
    /** connect to signal emitted before profiling */
    signal<void (uint64_t)>::connection on_prepend_profile(signal<void (uint64_t)>::slot_function_type const& slot);
    /** connect to signal emitted after profiling */
    signal<void (uint64_t)>::connection on_append_profile(signal<void (uint64_t)>::slot_function_type const& slot);
    /** profile */
    void profile(uint64_t step);
    /** Lua bindings */
    static void luaopen(lua_State* L);

private:
    /** pair type for accumulator with description */
    typedef std::pair<accumulator_type, std::string> accumulator_pair_type;

    /**
     * push copy of accumulator into container, and *reset* accumulator
     */
    void push(boost::shared_ptr<accumulator_type> acc, std::string const& desc);
    /**
     * write log entries for all runtime accumulators,
     * sorted by their total accumulated runtime
     */
    void log();

    /** signal emitted before profiling */
    signal<void (uint64_t)> on_prepend_profile_;
    /** signal emitted for profiling */
    signal_type on_profile_;
    /** signal emitted after profiling */
    signal<void (uint64_t)> on_append_profile_;
    /** accumulators with descriptions */
    std::vector<accumulator_pair_type> accumulators_;
};

} // namespace utility
} // namespace halmd

#endif /* ! HALMD_UTILITY_PROFILER_HPP */
