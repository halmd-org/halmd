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

#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <functional>

#include <halmd/io/logger.hpp>
#include <halmd/observables/sampler.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {

sampler::sampler(
    boost::shared_ptr<clock_type> clock
  , boost::shared_ptr<core_type> core
)
  : clock_(clock)
  , core_(core) {}

void sampler::setup()
{
    LOG("setting up simulation box");

    on_prepare_();
    core_->setup();
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

            {
                scoped_timer_type timer(runtime_.prepare);
                on_prepare_();
            }

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

/**
 * Forward signal to slot at given interval
 */
void sampler::prepare(std::function<void ()> const& slot, step_type interval) const
{
    if (clock_->step() % interval == 0) {
        slot();
    }
}

/**
 * Forward signal to slot at given interval
 */
void sampler::sample(std::function<void ()> const& slot, step_type interval) const
{
    if (clock_->step() % interval == 0) {
        slot();
    }
}

static std::function<void ()>
wrap_abort(boost::shared_ptr<mdsim::clock const> clock)
{
    return [=]() {
        throw std::runtime_error("gracefully aborting simulation at step " + boost::lexical_cast<std::string>(clock->step()));
    };
}

static std::function<void ()>
wrap_setup(boost::shared_ptr<sampler> self)
{
    return [=]() {
        return self->setup();
    };
}

static std::function<void (sampler::step_type)>
wrap_run(boost::shared_ptr<sampler> self)
{
    return [=](sampler::step_type steps) {
        return self->run(steps);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_start(boost::shared_ptr<sampler> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_start(slot);
    };
}

static std::function<void (std::function<void ()> const&)>
wrap_on_finish(boost::shared_ptr<sampler> self)
{
    return [=](std::function<void ()> const& slot) {
        return self->on_finish(slot);
    };
}

static std::function<void (std::function<void ()> const&, sampler::step_type)>
wrap_on_prepare(boost::shared_ptr<sampler> self)
{
    return [=](std::function<void ()> const& slot, sampler::step_type steps) {
        return self->on_prepare(slot, steps);
    };
}

static std::function<void (std::function<void ()> const&, sampler::step_type)>
wrap_on_sample(boost::shared_ptr<sampler> self)
{
    return [=](std::function<void ()> const& slot, sampler::step_type steps) {
        return self->on_sample(slot, steps);
    };
}

void sampler::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        class_<sampler, boost::shared_ptr<sampler> >("sampler")
            .def(constructor<
                boost::shared_ptr<sampler::clock_type>
              , boost::shared_ptr<sampler::core_type>
            >())
            .property("setup", &wrap_setup)
            .property("run", &wrap_run)
            .property("on_start", &wrap_on_start)
            .property("on_finish", &wrap_on_finish)
            .property("on_prepare", &wrap_on_prepare)
            .property("on_sample", &wrap_on_sample)
            .scope
            [
                class_<runtime>("runtime")
                    .def_readonly("total", &runtime::total)
                    .def_readonly("start", &runtime::start)
                    .def_readonly("prepare", &runtime::prepare)
                    .def_readonly("sample", &runtime::sample)
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
