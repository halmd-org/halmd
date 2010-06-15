/*
 * Copyright © 2008-2010  Peter Colberg
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
    modules::required<_Self, force_type>();
    modules::required<_Self, neighbour_type>();
    modules::optional<_Self, sort_type>();
    modules::required<_Self, integrator_type>();
    modules::required<_Self, position_type>();
    modules::required<_Self, velocity_type>();
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
core<dimension>::core(po::options const& vm)
  // dependency injection
  : force(modules::fetch<force_type>(vm))
  , neighbour(modules::fetch<neighbour_type>(vm))
  , sort(modules::fetch<sort_type>(vm))
  , integrator(modules::fetch<integrator_type>(vm))
  , position(modules::fetch<position_type>(vm))
  , velocity(modules::fetch<velocity_type>(vm))
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
    integrator->integrate();
    if (neighbour->check()) {
        if (sort) {
            sort->order();
        }
        neighbour->update();
    }
    force->compute();
    integrator->finalize();
}

// explicit instantiation
template class core<3>;
template class core<2>;

} // namespace mdsim

template class module<mdsim::core<3> >;
template class module<mdsim::core<2> >;

} // namespace halmd
