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

template <int dimension>
boost::shared_ptr<neighbor<dimension> >
factory<neighbor<dimension> >::fetch(options const& vm)
{
    if (!neighbor_) {
#ifdef USE_HOST_SINGLE_PRECISION
        neighbor_.reset(new host::neighbor<dimension, float>(vm));
#else
        neighbor_.reset(new host::neighbor<dimension, double>(vm));
#endif
    }
    return neighbor_;
}

template <> factory<neighbor<3> >::neighbor_ptr factory<neighbor<3> >::neighbor_ = neighbor_ptr();
template class factory<neighbor<3> >;
template <> factory<neighbor<2> >::neighbor_ptr factory<neighbor<2> >::neighbor_ = neighbor_ptr();
template class factory<neighbor<2> >;

}} // namespace halmd::mdsim
