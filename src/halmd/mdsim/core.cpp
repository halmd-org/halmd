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

#include <cmath>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/util/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Assemble module options
 */
template <int dimension>
void core<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("dimension", po::value<int>()->default_value(3),
         "dimension of positional coordinates")
        ;
}

/**
 * Resolve module dependencies
 */
template <int dimension>
void core<dimension>::depends()
{
    modules::depends<_Self, force_type>::required();
    modules::depends<_Self, neighbour_type>::required();
    modules::depends<_Self, sort_type>::optional();
    modules::depends<_Self, integrator_type>::required();
    modules::depends<_Self, position_type>::required();
    modules::depends<_Self, velocity_type>::required();
}

template <int dimension>
void core<dimension>::select(po::options const& vm)
{
    if (vm["dimension"].as<int>() != dimension) {
        throw unsuitable_module("mismatching option dimension");
    }
}

/**
 * Initialize simulation
 */
template <int dimension>
core<dimension>::core(modules::factory& factory, po::options const& vm)
  // dependency injection
  : force(modules::fetch<force_type>(factory, vm))
  , neighbour(modules::fetch<neighbour_type>(factory, vm))
  , sort(modules::fetch<sort_type>(factory, vm))
  , integrator(modules::fetch<integrator_type>(factory, vm))
  , position(modules::fetch<position_type>(factory, vm))
  , velocity(modules::fetch<velocity_type>(factory, vm))
{
    position->set();
    velocity->set();
    neighbour->update();
    force->compute();
}

/**
 * Perform a single MD integration step
 */
template <int dimension>
inline void core<dimension>::mdstep()
{
    boost::array<high_resolution_timer, 8> timer;
    timer[0].record();

    integrator->integrate();
    timer[1].record();

    bool nbr_check = neighbour->check();
    timer[2].record();
    if (nbr_check) {
        if (sort) {
            sort->order();
        }
        timer[3].record();
        neighbour->update();
    }
    timer[4].record();

    force->compute();
    timer[5].record();

    integrator->finalize();
    timer[6].record();

    runtimes["mdstep"] += timer[6] - timer[0];
    runtimes["integration"] += (timer[1] - timer[0]) + (timer[6] - timer[5]);
    runtimes["update_forces"] += timer[5] - timer[4];
    runtimes["check_neighbour_update"] += timer[2] - timer[1];
    if (nbr_check) {
        if (sort) {
            runtimes["particle_sort"] += timer[3] - timer[2];
        }
        runtimes["update_neighbours"] += timer[4] - timer[3];
    }
    timer[7].record();

    runtimes["update_timers"] += timer[7] - timer[6];
}

// explicit instantiation
template class core<3>;
template class core<2>;

} // namespace mdsim

template class module<mdsim::core<3> >;
template class module<mdsim::core<2> >;

} // namespace halmd
