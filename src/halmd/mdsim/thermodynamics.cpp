/*
 * Copyright © 2010  Felix Höfling
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

#include <halmd/mdsim/thermodynamics.hpp>
#include <halmd/utility/module.hpp>

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
void thermodynamics<dimension>::depends()
{
    modules::depends<_Self, box_type>::required();
}

template <int dimension>
thermodynamics<dimension>::thermodynamics(modules::factory& factory, po::options const& vm)
  // dependency injection
  : box(modules::fetch<box_type>(factory, vm))
{
}

// explicit instantiation
template class thermodynamics<3>;
template class thermodynamics<2>;

} // namespace mdsim

} // namespace halmd
