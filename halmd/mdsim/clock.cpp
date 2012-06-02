/*
 * Copyright © 2011-2012 Peter Colberg
 * Copyright © 2011 Felix Höfling
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

namespace halmd {
namespace mdsim {

clock::clock() : step_(0), time_(0), step_origin_(0), time_origin_(0) {}

clock::time_type clock::timestep() const
{
    if (!timestep_) {
        throw std::logic_error("time step has not been set");
    }
    return *timestep_;
}

void clock::set_timestep(time_type timestep)
{
    step_origin_ = step_;
    time_origin_ = time_;
    timestep_ = timestep;

    // propagate timestep to integrator(s)
    on_set_timestep_(*timestep_);

    LOG("integration time step: " << *timestep_);
}

static std::function<void (clock::time_type)>
wrap_set_timestep(boost::shared_ptr<clock> self)
{
    return [=](clock::time_type timestep) {
        self->set_timestep(timestep);
    };
}

static std::function<connection (std::function<void (clock::time_type)> const&)>
wrap_on_set_timestep(boost::shared_ptr<clock> self)
{
    return [=](std::function<void (clock::time_type)> const& slot) {
        return self->on_set_timestep(slot);
    };
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_clock(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<clock, boost::shared_ptr<clock> >("clock")
                .def(constructor<>())
                .property("set_timestep", &wrap_set_timestep)
                .property("on_set_timestep", &wrap_on_set_timestep)
                .property("step", &clock::step)
                .property("time", &clock::time)
                .property("timestep", &clock::timestep)
        ]
    ];
    return 0;
}

} // namespace mdsim
} // namespace halmd
