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
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables
{

/**
 * Initialize simulation
 */
template <int dimension>
sampler<dimension>::sampler(
    shared_ptr<core_type> core
  , uint64_t steps
)
  : core_(core)
  , steps_(steps)
  , total_time_(steps_ * core_->integrator->timestep())
{
    LOG("number of integration steps: " << steps_);
    LOG("integration time: " << total_time_);
}

/**
 * register module runtime accumulators
 */
template <int dimension>
void sampler<dimension>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.total, "total", "total simulation runtime");
}

/**
 * Run simulation
 */
template <int dimension>
void sampler<dimension>::run()
{

    LOG("setting up simulation box");

    on_prepare_(core_->time());
    core_->prepare();
    on_sample_(core_->time());

    on_start_(core_->time());

    LOG("starting simulation run");
    {
        scoped_timer<timer> timer_(runtime_.total);

        while (core_->step_counter() < steps_) {
            on_prepare_(core_->time());

            // perform complete MD integration step
            core_->mdstep();

            on_sample_(core_->time());
        }
    }
    LOG("finished simulation run");

    on_finish_(core_->time());
}

/**
 * Connect slot to signal emitted before starting simulation run
 */
template <int dimension>
void sampler<dimension>::on_start(slot_function_type const& slot)
{
    on_start_.connect(slot);
}

/**
 * Connect slot to signal emitted before MD integration step
 */
template <int dimension>
void sampler<dimension>::on_prepare(slot_function_type const& slot, uint64_t interval)
{
    on_prepare_.connect(
        bind(&sampler::prepare, this, slot, interval, _1)
    );
}

/**
 * Connect slot to signal emitted after MD integration step
 */
template <int dimension>
void sampler<dimension>::on_sample(slot_function_type const& slot, uint64_t interval)
{
    on_sample_.connect(
        bind(&sampler::sample, this, slot, interval, _1)
    );
}

/**
 * Connect slot to signal emitted after finishing simulation run
 */
template <int dimension>
void sampler<dimension>::on_finish(slot_function_type const& slot)
{
    on_finish_.connect(slot);
}

/**
 * Forward signal to slot at given interval
 */
template <int dimension>
void sampler<dimension>::prepare(slot_function_type const& slot, uint64_t interval, double time) const
{
    uint64_t step = core_->step_counter(); // FIXME +1 if not core->prepare()
    if (step == 0 || step % interval == 0 || step == steps_) {
        slot(time);
    }
}

/**
 * Forward signal to slot at given interval
 */
template <int dimension>
void sampler<dimension>::sample(slot_function_type const& slot, uint64_t interval, double time) const
{
    uint64_t step = core_->step_counter();
    if (step == 0 || step % interval == 0 || step == steps_) {
        slot(time);
    }
}

template <int dimension>
void sampler<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("sampler_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            class_<sampler, shared_ptr<runner>, runner>(class_name.c_str())
                .def(constructor<
                    shared_ptr<core_type>
                  , uint64_t
                >())
                .def("register_runtimes", &sampler::register_runtimes)
                .property("steps", &sampler::steps)
                .property("total_time", &sampler::total_time)
                .def("on_start", &sampler::on_start)
                .def("on_finish", &sampler::on_finish)
                .def("on_prepare", &sampler::on_prepare)
                .def("on_sample", &sampler::on_sample)
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &sampler<3>::luaopen
    ]
    [
        &sampler<2>::luaopen
    ];
}

} // namespace

// explicit instantiation
template class sampler<3>;
template class sampler<2>;

} // namespace observables

} // namespace halmd
