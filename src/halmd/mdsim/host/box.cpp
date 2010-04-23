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

#include <algorithm>
#include <cmath>
#include <numeric>

#include <halmd/mdsim/host/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/util/logger.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host
{

/**
 * Set box edge lengths
 */
template <int dimension>
box<dimension>::box(options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(module<particle_type>::fetch(vm))
  // initialize parameters
  , length_half_(.5 * length_)
{}

template <int dimension>
typename box<dimension>::pointer
box<dimension>::create(options const& vm)
{
#ifdef USE_HOST_SINGLE_PRECISION
    if (module<host::particle<dimension, float> >::fetch(vm)) {
        return pointer(new box<dimension>(vm));
    }
#else
    if (module<host::particle<dimension, double> >::fetch(vm)) {
        return pointer(new box<dimension>(vm));
    }
#endif
    return pointer();
}

// explicit instantiation
template class box<3>;
template class box<2>;

}} // namespace mdsim::host

template class module<mdsim::host::box<3> >;
template class module<mdsim::host::box<2> >;

} // namespace halmd
