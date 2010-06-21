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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/integrator.hpp>

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
void integrator<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("integrator", po::value<string>()->default_value("verlet"),
         "specify integration module")
        ;

    po::options_description group("Integrator");
    group.add_options()
        ("timestep,h", po::value<double>()->default_value(0.001),
         "integration timestep")
        ;
    desc.add(group);
}

template <int dimension>
integrator<dimension>::integrator(modules::factory& factory, po::options const& vm)
  // set parameters
  : timestep_(vm["timestep"].as<double>())
{
    LOG("integration timestep: " << timestep_);
}

// explicit instantiation
template class integrator<3>;
template class integrator<2>;

} // namespace mdsim

} // namespace halmd
