/*
 * Copyright © 2010-2011  Felix Höfling and Peter Colberg
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

#include <halmd/io/logger.hpp>
#include <halmd/observables/sampler.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {

/**
 * Initialize simulation
 */
sampler::sampler(
    shared_ptr<clock_type const> clock
  , shared_ptr<core_type> core
  , step_type steps
)
  : clock_(clock)
  , core_(core)
  , steps_(steps)
  , total_time_(steps_ * clock_->timestep())
{
    LOG("number of integration steps: " << steps_);
    LOG("integration time: " << total_time_);
}

/**
 * register module runtime accumulators
 */
void sampler::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.total, "total", "total simulation runtime");
}

/**
 * Run simulation
 */
void sampler::run()
{
    LOG("setting up simulation box");

    on_prepare_(clock_->step());
    core_->setup();
    on_sample_(clock_->step());

    on_start_(clock_->step());

    LOG("starting simulation run");
    {
        scoped_timer<timer> timer_(runtime_.total);

        while (clock_->step() < steps_) {
            on_prepare_(clock_->step() + 1); //< step counter is increased by call to mdstep()

            // perform complete MD integration step
            core_->mdstep();

            on_sample_(clock_->step());
        }
    }
    LOG("finished simulation run");

    on_finish_(clock_->step());
}

/**
 * Connect slot to signal emitted before starting simulation run
 */
sampler::connection_type sampler::on_start(slot_function_type const& slot)
{
    return on_start_.connect(slot);
}

/**
 * Connect slot to signal emitted before MD integration step
 */
sampler::connection_type sampler::on_prepare(slot_function_type const& slot, step_type interval)
{
    return on_prepare_.connect(bind(&sampler::prepare, this, slot, interval, _1));
}

/**
 * Connect slot to signal emitted after MD integration step
 */
sampler::connection_type sampler::on_sample(slot_function_type const& slot, step_type interval)
{
    return on_sample_.connect(bind(&sampler::sample, this, slot, interval, _1));
}

/**
 * Connect slot to signal emitted after finishing simulation run
 */
sampler::connection_type sampler::on_finish(slot_function_type const& slot)
{
    return on_finish_.connect(slot);
}

/**
 * Forward signal to slot at given interval
 */
void sampler::prepare(slot_function_type const& slot, step_type interval, step_type step) const
{
    if (step == 0 || step % interval == 0 || step == steps_) {
        slot(step);
    }
}

/**
 * Forward signal to slot at given interval
 */
void sampler::sample(slot_function_type const& slot, step_type interval, step_type step) const
{
    if (step == 0 || step % interval == 0 || step == steps_) {
        slot(step);
    }
}

HALMD_LUA_API int luaopen_libhalmd_observables_sampler(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        class_<sampler, shared_ptr<sampler> >("sampler")
            .def(constructor<
                shared_ptr<sampler::clock_type const>
              , shared_ptr<sampler::core_type>
              , sampler::step_type
            >())
            .def("register_runtimes", &sampler::register_runtimes)
            .def("on_start", &sampler::on_start)
            .def("on_finish", &sampler::on_finish)
            .def("on_prepare", &sampler::on_prepare)
            .def("on_sample", &sampler::on_sample)
            .def("run", &sampler::run)
            .property("steps", &sampler::steps)
            .property("total_time", &sampler::total_time)
    ];
    return 0;
}

} // namespace observables
} // namespace halmd
