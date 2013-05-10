/*
 * Copyright © 2010-2012 Peter Colberg
 * Copyright © 2010-2011 Felix Höfling
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

#include <exception>
#include <functional>

#include <halmd/io/logger.hpp>
#include <halmd/observables/sampler.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {

sampler::sampler(
    std::shared_ptr<clock_type> clock
  , std::shared_ptr<core_type> core
)
  : clock_(clock)
  , core_(core) {}

void sampler::sample()
{
    LOG_TRACE("sample state at step " << clock_->step());

    on_sample_();
}

void sampler::run(step_type steps)
{
    {
        scoped_timer_type timer(runtime_.start);
        on_start_();
    }

    LOG("starting simulation run over " << steps << " integration steps");
    {
        scoped_timer_type timer(runtime_.total);

        step_type limit = clock_->step() + steps;

        while (clock_->step() < limit) {
            // increment 1-based simulation step
            clock_->advance();

            LOG_TRACE("performing MD step #" << clock_->step());

            // perform complete MD integration step
            core_->mdstep();

            {
                scoped_timer_type timer(runtime_.sample);
                on_sample_();
            }
        }
    }
    LOG("finished simulation run");

    {
        scoped_timer_type timer(runtime_.finish);
        on_finish_();
    }
}

connection sampler::on_sample(std::function<void ()> const& slot, step_type interval, step_type start)
{
    if (interval == 0) {
        throw std::logic_error("Slot must not be connected to signal 'on_sample' with zero sampling interval");
    }
    return on_sample_.connect([=]() {
        step_type step = clock_->step();
        if(step >= start && (step - start) % interval == 0) {
            slot();
        }
    });
}

connection sampler::on_start(std::function<void ()> const& slot)
{
    return on_start_.connect(slot);
}

connection sampler::on_finish(std::function<void ()> const& slot)
{
    return on_finish_.connect(slot);
}

static std::function<void ()>
wrap_abort(std::shared_ptr<mdsim::clock const> clock)
{
    return [=]() {
        throw std::runtime_error("gracefully aborting simulation at step " + std::to_string(clock->step()));
    };
}

void sampler::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        class_<sampler, std::shared_ptr<sampler> >("sampler")
            .def(constructor<
                std::shared_ptr<sampler::clock_type>
              , std::shared_ptr<sampler::core_type>
            >())
            .def("sample", &sampler::sample)
            .def("run", &sampler::run)
            .def("on_sample", &sampler::on_sample)
            .def("on_start", &sampler::on_start)
            .def("on_finish", &sampler::on_finish)
            .scope
            [
                class_<runtime>("runtime")
                    .def_readonly("total", &runtime::total)
                    .def_readonly("sample", &runtime::sample)
                    .def_readonly("start", &runtime::start)
                    .def_readonly("finish", &runtime::finish)

              , def("abort", &wrap_abort)
            ]
            .def_readonly("runtime", &sampler::runtime_)
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_sampler(lua_State* L)
{
    sampler::luaopen(L);
    return 0;
}

} // namespace observables
} // namespace halmd
