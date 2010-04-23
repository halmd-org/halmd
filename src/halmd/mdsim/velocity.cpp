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

#include <halmd/mdsim/host/velocity/boltzmann.hpp>
#include <halmd/mdsim/velocity.hpp>
#include <halmd/util/logger.hpp>

namespace halmd
{
namespace mdsim
{

template <int dimension>
typename velocity<dimension>::pointer
velocity<dimension>::create(options const& vm)
{
#ifdef USE_HOST_SINGLE_PRECISION
    return pointer(new host::velocity::boltzmann<dimension, float>(vm));
#else
    return pointer(new host::velocity::boltzmann<dimension, double>(vm));
#endif
}

template class velocity<3>;
template class velocity<2>;

} // namespace mdsim

template class module<mdsim::velocity<3> >;
template class module<mdsim::velocity<2> >;

} // namespace halmd
