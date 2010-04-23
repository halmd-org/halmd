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

#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/host/forces/lj.hpp>
#include <halmd/util/logger.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

template <int dimension>
force<dimension>::force(options const& vm)
  // dependency injection
  : particle(module<particle_type>::fetch(vm))
  // allocate result variables
  , virial_(particle->ntype)
{
}

template <int dimension>
typename force<dimension>::pointer
force<dimension>::create(options const& vm)
{
#ifdef USE_HOST_SINGLE_PRECISION
    return pointer(new host::forces::lj<dimension, float>(vm));
#else
    return pointer(new host::forces::lj<dimension, double>(vm));
#endif
}

// explicit instantiation
template class force<3>;
template class force<2>;

} // namespace mdsim

template class module<mdsim::force<3> >;
template class module<mdsim::force<2> >;

} // namespace halmd
