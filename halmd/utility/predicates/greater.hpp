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

#ifndef HALMD_UTILITY_PREDICATES_GREATER_HPP
#define HALMD_UTILITY_PREDICATES_GREATER_HPP

#include <boost/function.hpp>
#include <lua.hpp>

#include <halmd/utility/signal.hpp>

namespace halmd {
namespace predicates {

template <typename value_type>
class greater
{
public:
    typedef halmd::signal<void ()> signal_type;
    typedef signal_type::slot_function_type slot_function_type;
    typedef boost::function<value_type ()> function_type;

    greater(function_type const& func, value_type const& value) : func_(func), value_(value) {}

    connection on_greater(slot_function_type const& slot)
    {
        return on_greater_.connect(slot);
    }

    void evaluate() const
    {
        if (func_() > value_) {
            on_greater_();
        }
    }

    static void luaopen(lua_State* L, char const* class_name);

private:
    function_type func_;
    value_type value_;
    signal_type on_greater_;
};

} // namespace predicates
} // namespace halmd

#endif /* ! HALMD_UTILITY_PREDICATES_GREATER_HPP */
