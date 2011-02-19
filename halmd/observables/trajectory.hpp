/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_TRAJECTORY_HPP
#define HALMD_OBSERVABLES_TRAJECTORY_HPP

#include <boost/function.hpp>
#include <lua.hpp>

namespace halmd
{
namespace observables
{

template <int dimension>
class trajectory
{
public:
    static void luaopen(lua_State* L);

    trajectory() {}
    virtual ~trajectory() {};
    virtual void acquire(double time) = 0;

    // data stream interface
    virtual void register_request(uint64_t step, boost::function<void(uint64_t)> callback) = 0;
    virtual void notify(uint64_t step) = 0;
};

} // namespace observables

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_TRAJECTORY_HPP */
