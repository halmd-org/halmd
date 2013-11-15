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

#ifndef HALMD_MDSIM_CORE_HPP
#define HALMD_MDSIM_CORE_HPP

#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

/** HAL’s MD package */
namespace halmd {
/** Molecular Dynamics simulation modules */
namespace mdsim {

class core
{
public:
    typedef halmd::signal<void ()> signal_type;
    typedef signal_type::slot_function_type slot_function_type;

    void mdstep();

    connection on_prepend_integrate(slot_function_type const& slot)
    {
        return on_prepend_integrate_.connect(slot);
    }

    connection on_integrate(slot_function_type const& slot)
    {
        return on_integrate_.connect(slot);
    }

    connection on_append_integrate(slot_function_type const& slot)
    {
        return on_append_integrate_.connect(slot);
    }

    connection on_prepend_finalize(slot_function_type const& slot)
    {
        return on_prepend_finalize_.connect(slot);
    }

    connection on_finalize(slot_function_type const& slot)
    {
        return on_finalize_.connect(slot);
    }

    connection on_append_finalize(slot_function_type const& slot)
    {
        return on_append_finalize_.connect(slot);
    }

    /** Lua bindings */
    static void luaopen(lua_State* L);

private:
    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type mdstep;
    };

    signal_type on_prepend_integrate_;
    signal_type on_integrate_;
    signal_type on_append_integrate_;
    signal_type on_prepend_finalize_;
    signal_type on_finalize_;
    signal_type on_append_finalize_;
    runtime runtime_;
};

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_CORE_HPP */
