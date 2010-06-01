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

#include <halmd/mdsim/force.hpp>

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
void force<dimension>::resolve(po::options const& vm)
{
    module<particle_type>::required(vm);
}

template <int dimension>
force<dimension>::force(po::options const& vm)
  // dependency injection
  : particle(module<particle_type>::fetch(vm))
{
}

// explicit instantiation
template class force<3>;
template class force<2>;

} // namespace mdsim

} // namespace halmd
