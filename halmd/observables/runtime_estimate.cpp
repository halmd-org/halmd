/*
 * Copyright © 2008-2011  Felix Höfling and Peter Colberg
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
#include <iomanip>
#include <sstream>

#include <halmd/io/logger.hpp>
#include <halmd/observables/runtime_estimate.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {

runtime_estimate::runtime_estimate(
    boost::shared_ptr<clock_type> clock
  , step_type total_steps
)
  // dependency injection
  : clock_(clock)
  // initialise members
  , step_start_(clock_->step())
  , step_stop_(step_start_ + total_steps)
  // A circular buffer of two timer/step samples allows estimating the
  // runtime in between calls to runtime_estimate::sample, e.g. triggered
  // by the user through a POSIX signal. For the estimate, the front of
  // the circular buffer will be used, thus the measured time interval
  // is at least as long as the sample interval.
  , timer_(2)
{
}

void runtime_estimate::sample()
{
    timer_.push_back(make_pair(timer(), clock_->step()));
}

void runtime_estimate::estimate() const
{
    // linear extrapolation of measured steps to total number of steps
    timer_pair_type const& timer = timer_.front();
    step_type const& step = clock_->step();
    double eta = (timer.first.elapsed() / (step - timer.second)) * (step_stop_ - step);
    LOG(format_time(eta, 1) << " estimated remaining runtime at step " << step);
}

string runtime_estimate::format_time(double time, unsigned int prec)
{
    ostringstream os;
    if (time < 60)
        os << std::fixed << std::setprecision(prec) << time << " s";
    else if (time < 3600)
        os << std::fixed << std::setprecision(prec) << (time / 60) << " min";
    else if (time < 86400)
        os << std::fixed << std::setprecision(prec) << (time / 3600) << " h";
    else
        os << std::fixed << std::setprecision(prec) << (time / 86400) << " d";
    return os.str();
}

static runtime_estimate::slot_function_type
wrap_sample(boost::shared_ptr<runtime_estimate> instance)
{
    return bind(&runtime_estimate::sample, instance);
}

static runtime_estimate::slot_function_type
wrap_estimate(boost::shared_ptr<runtime_estimate> instance)
{
    return bind(&runtime_estimate::estimate, instance);
}

void runtime_estimate::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<runtime_estimate, boost::shared_ptr<runtime_estimate> >("runtime_estimate")
                .def(constructor<
                    boost::shared_ptr<clock_type>
                  , step_type
                >())
                .property("sample", &wrap_sample)
                .property("estimate", &wrap_estimate)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_runtime_estimate(lua_State* L)
{
    runtime_estimate::luaopen(L);
    return 0;
}

} // namespace observables
} // namespace halmd
