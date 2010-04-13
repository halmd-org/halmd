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

#include <halmd/mdsim/host/neighbor.hpp>
#include <halmd/mdsim/neighbor.hpp>
#include <halmd/util/log.hpp>

namespace halmd { namespace mdsim
{

template <int dimension, typename float_type>
boost::shared_ptr<neighbor<dimension, float_type> >
factory<neighbor<dimension, float_type> >::fetch(options const& vm)
{
    if (!neighbor_) {
        neighbor_.reset(new host::neighbor<dimension, float_type>(vm));
    }
    return neighbor_;
}

#ifndef USE_HOST_SINGLE_PRECISION
template <> factory<neighbor<3, double> >::neighbor_ptr factory<neighbor<3, double> >::neighbor_ = neighbor_ptr();
template class factory<neighbor<3, double> >;
template <> factory<neighbor<2, double> >::neighbor_ptr factory<neighbor<2, double> >::neighbor_ = neighbor_ptr();
template class factory<neighbor<2, double> >;
#else
template <> factory<neighbor<3, float> >::neighbor_ptr factory<neighbor<3, float> >::neighbor_ = neighbor_ptr();
template class factory<neighbor<3, float> >;
template <> factory<neighbor<2, float> >::neighbor_ptr factory<neighbor<2, float> >::neighbor_ = neighbor_ptr();
template class factory<neighbor<2, float> >;
#endif

}} // namespace halmd::mdsim
