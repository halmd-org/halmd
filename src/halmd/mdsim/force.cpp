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

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim
{

template <int dimension, typename float_type>
force<dimension, float_type>::force(particle_ptr particle)
    // dependency injection
    : particle(dynamic_pointer_cast<particle_type>(particle))
    // allocate result variables
    , virial_(particle->ntype)
{
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class force<3, double>;
template class force<2, double>;
#else
template class force<3, float>;
template class force<2, float>;
#endif

}} // namespace halmd::mdsim
