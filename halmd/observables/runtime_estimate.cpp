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

runtime_estimate::runtime_estimate(shared_ptr<clock_type> clock, uint64_t total_steps, uint64_t step_start)
  // dependency injection
  : clock(clock)
  // initialise members
  , step_start_(step_start)
  , step_stop_(step_start + total_steps)
{
    gettimeofday(&start_time_, NULL);
}

/**
 * trigger estimate of remaining runtime and output to log file
 */
void runtime_estimate::sample() const
{
    uint64_t step = clock->step();
    if (step > step_start_) {
        double eta = value(step);
        LOG(format_time(eta, 1) << " estimated remaining runtime at step " << step);
    }
}

/**
 * estimate remaining runtime
 */
double runtime_estimate::value(uint64_t step) const
{
    // return if called at initial step
    if (step == step_start_) {
        return 0;
    }

    // determine passed wall clock time since start
    timeval stop_time;
    gettimeofday(&stop_time, NULL);

    timeval tv;
    timersub(&stop_time, &start_time_, &tv);

    // linear extrapolation
    return (tv.tv_sec + tv.tv_usec * 1e-6) * (step_stop_ - step) / (step - step_start_);
}

/**
 * format time given in seconds
 */
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


runtime_estimate::slot_function_type
sample_wrapper(shared_ptr<runtime_estimate> instance)
{
    return bind(&runtime_estimate::sample, instance);
}

void runtime_estimate::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("runtime_estimate_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<runtime_estimate, shared_ptr<runtime_estimate> >(class_name.c_str())
                .def(constructor<
                    shared_ptr<clock_type>
                  , uint64_t, uint64_t
                >())
                .property("sample", &sample_wrapper)
                .property("value", &runtime_estimate::value)
                .def("on_sample", &runtime_estimate::on_sample)
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
