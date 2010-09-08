/*
 * Copyright © 2010  Peter Colberg and Felix Höfling
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

#include <algorithm>
#include <boost/bind.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/scripts/thermostat.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace scripts
{

/**
 * Assemble module options
 */
template <int dimension>
void thermostat<dimension>::options(po::options_description& desc)
{
    po::options_description group("Simulation");
    group.add_options()
        ("thermostat", po::value<float>()->required(),
         "heat bath collision rate")
        ;
    desc.add(group);
}

/**
 * Resolve module dependencies
 */
template <int dimension>
void thermostat<dimension>::depends()
{
    // FIXME ugly hack to override base module
    modules::depends<_Self, _Self>::required();
}

template <int dimension>
thermostat<dimension>::thermostat(modules::factory& factory, po::options const& vm)
  : _Base(factory, vm)
  // set parameters
  , rate_(vm["thermostat"].as<float>())
  , period_(max(static_cast<int>(round(1 / (rate_ * core->integrator->timestep()))), 1))
{
    LOG("heat bath collision rate: " << rate_);
    LOG("heat bath coupling period: " << period_);
}

/**
 * Run simulation
 */
template <int dimension>
void thermostat<dimension>::run()
{
    core->prepare();
    sampler->sample();

    LOG("starting thermostat run");

    while (core->step_counter() < _Base::steps_) {
        // perform complete MD integration step
        core->mdstep();

        // assign new random velocities (Andersen thermostat)
        if (!(core->step_counter() % period_)) {
            core->velocity->set();
        }

        // sample system state and properties
        sampler->sample();
    }

    LOG("finished thermostat run");

    for_each(
        profile_writers.begin()
      , profile_writers.end()
      , bind(&profile_writer_type::write, _1)
    );
}

// explicit instantiation
template class thermostat<3>;
template class thermostat<2>;

} // namespace scripts

template class module<scripts::thermostat<3> >;
template class module<scripts::thermostat<2> >;

} // namespace halmd
