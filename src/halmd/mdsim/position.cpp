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

#include <halmd/mdsim/host/position/lattice.hpp>
#include <halmd/mdsim/position.hpp>
#include <halmd/util/log.hpp>

namespace halmd { namespace mdsim
{

template <int dimension>
boost::shared_ptr<position<dimension> >
factory<position<dimension> >::fetch(options const& vm)
{
    if (!position_) {
#ifdef USE_HOST_SINGLE_PRECISION
        position_.reset(new host::position::lattice<dimension, float>(vm));
#else
        position_.reset(new host::position::lattice<dimension, double>(vm));
#endif
    }
    return position_;
}

template <> factory<position<3> >::position_ptr factory<position<3> >::position_ = position_ptr();
template class factory<position<3> >;
template <> factory<position<2> >::position_ptr factory<position<2> >::position_ = position_ptr();
template class factory<position<2> >;

}} // namespace halmd::mdsim
