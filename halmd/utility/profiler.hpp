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

#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd {
namespace utility {

/**
 * Module runtime profiler
 *
 * This class manages the runtime accumulators scattered over the
 * HALMD modules. It allows granular profiling, i.e. only a subset
 * of accumulators may be registered with the profiler.
 *
 * In the modules, accumulators are defined as member variables in a
 * runtime struct, and exported with Luabind's .def_readonly(). The
 * profiler Lua module then calls the profile() function of a module,
 * which calls profiler:on_profile for each accumulator along with
 * a description.
 *
 * Profiler holds the accumulators as shared pointers. Lua automatically
 * creates a shared pointer of an accumulator when calling on_profile,
 * and Lua will keep alive the module which contains the accumulator
 * as long as the shared pointer is held by profiler.
 *
 * We allow disconnection of accumulators by use of the halmd::slots
 * container, which returns a connection object when an accumulator
 * is connected (i.e. inserted).
 */
class profiler
{
private:
    typedef signal<void ()> signal_type;

public:
    typedef timer timer_type;
    typedef accumulator<double> accumulator_type;
    typedef scoped_timer<timer_type> scoped_timer_type;
    typedef signal_type::slot_function_type slot_function_type;

    /** logs timer resolution */
    profiler();
    /** connect accumulator for profiling */
    connection on_profile(std::shared_ptr<accumulator_type> acc, std::string const& desc);
    /** connect to signal emitted before profiling */
    connection on_prepend_profile(slot_function_type const& slot);
    /** connect to signal emitted after profiling */
    connection on_append_profile(slot_function_type const& slot);
    /** log and reset runtime accumulators */
    void profile();
    /** Lua bindings */
    static void luaopen(lua_State* L);

private:
    /** accumulator with description */
    typedef std::pair<std::shared_ptr<accumulator_type>, std::string> accumulator_pair_type;
    typedef slots<accumulator_pair_type> slots_type;
    typedef slots_type::const_iterator slots_const_iterator;

    void log() const;

    /** accumulators slots */
    slots_type accumulators_;
    /** signal emitted before profiling */
    signal_type on_prepend_profile_;
    /** signal emitted after profiling */
    signal_type on_append_profile_;
};

} // namespace utility
} // namespace halmd

#endif /* ! HALMD_UTILITY_PROFILER_HPP */
