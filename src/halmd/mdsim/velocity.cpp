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
#include <halmd/util/log.hpp>

namespace halmd { namespace mdsim
{

template <int dimension>
boost::shared_ptr<velocity<dimension> >
factory<velocity<dimension> >::fetch(options const& vm)
{
    if (!velocity_) {
#ifdef USE_HOST_SINGLE_PRECISION
        velocity_.reset(new host::velocity::boltzmann<dimension, float>(vm));
#else
        velocity_.reset(new host::velocity::boltzmann<dimension, double>(vm));
#endif
    }
    return velocity_;
}

template <> factory<velocity<3> >::velocity_ptr factory<velocity<3> >::velocity_ = velocity_ptr();
template class factory<velocity<3> >;
template <> factory<velocity<2> >::velocity_ptr factory<velocity<2> >::velocity_ = velocity_ptr();
template class factory<velocity<2> >;

}} // namespace halmd::mdsim
