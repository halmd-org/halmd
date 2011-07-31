/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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

#include <boost/shared_ptr.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {

clock::clock(time_type timestep)
  // initialise attributes
  : step_(0)
  , time_(0)
  , timestep_(timestep)
{
    LOG("simulation time-step: " << timestep_);
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_clock(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<clock, shared_ptr<clock> >("clock")
                .def(constructor<clock::time_type>())
                .property("step", &clock::step)
                .property("time", &clock::time)
                .property("timestep", &clock::timestep)
        ]
    ];
    return 0;
}

} // namespace mdsim
} // namespace halmd
