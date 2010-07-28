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
#include <halmd/mdsim/position.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Resolve module dependencies
 */
template <int dimension>
void position<dimension>::depends()
{
    modules::depends<_Self, profiler_type>::required();
}

/**
 * Assemble module options
 */
template <int dimension>
void position<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("position",
         po::value<string>()->default_value("lattice"),
         "initial particle positions module")
        ;
}

template <int dimension>
position<dimension>::position(modules::factory& factory, po::options const& vm)
  // dependency injection
  : profiler(modules::fetch<profiler_type>(factory, vm))
{}

template class position<3>;
template class position<2>;

} // namespace mdsim

} // namespace halmd
