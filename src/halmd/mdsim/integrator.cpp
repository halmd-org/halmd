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
#include <halmd/util/logger.hpp>

namespace halmd { namespace mdsim
{

template <int dimension>
typename module<integrator<dimension> >::pointer
module<integrator<dimension> >::fetch(options const& vm)
{
    if (!singleton_) {
#ifdef USE_HOST_SINGLE_PRECISION
        singleton_.reset(new host::integrator::verlet<dimension, float>(vm));
#else
        singleton_.reset(new host::integrator::verlet<dimension, double>(vm));
#endif
    }
    return singleton_;
}

template <int dimension> typename module<integrator<dimension> >::pointer module<integrator<dimension> >::singleton_;

template class module<integrator<3> >;
template class module<integrator<2> >;

}} // namespace halmd::mdsim
