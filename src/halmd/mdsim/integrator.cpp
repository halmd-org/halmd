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

#include <halmd/mdsim/host/integrator/verlet.hpp>
#include <halmd/mdsim/integrator.hpp>
#include <halmd/util/log.hpp>

namespace halmd { namespace mdsim
{

template <int dimension, typename float_type>
boost::shared_ptr<integrator<dimension, float_type> >
factory<integrator<dimension, float_type> >::fetch(options const& vm)
{
    if (!integrator_) {
        integrator_.reset(new host::integrator::verlet<dimension, float_type>(vm));
    }
    return integrator_;
}

#ifndef USE_HOST_SINGLE_PRECISION
template <> factory<integrator<3, double> >::integrator_ptr factory<integrator<3, double> >::integrator_ = integrator_ptr();
template class factory<integrator<3, double> >;
template <> factory<integrator<2, double> >::integrator_ptr factory<integrator<2, double> >::integrator_ = integrator_ptr();
template class factory<integrator<2, double> >;
#else
template <> factory<integrator<3, float> >::integrator_ptr factory<integrator<3, float> >::integrator_ = integrator_ptr();
template class factory<integrator<3, float> >;
template <> factory<integrator<2, float> >::integrator_ptr factory<integrator<2, float> >::integrator_ = integrator_ptr();
template class factory<integrator<2, float> >;
#endif

}} // namespace halmd::mdsim
