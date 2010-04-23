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

#include <halmd/mdsim/gpu/sample/trajectory.hpp>
#include <halmd/mdsim/samples/gpu/trajectory.hpp>
#include <halmd/util/logger.hpp>

namespace halmd
{
namespace mdsim { namespace samples { namespace gpu
{

template <int dimension, typename float_type>
typename trajectory<dimension, float_type>::pointer
trajectory<dimension, float_type>::create(options const& vm)
{
    return pointer(new mdsim::gpu::sample::trajectory<trajectory<dimension, float_type> >(vm));
}

template class trajectory<3, double>;
template class trajectory<2, double>;

}}} // namespace mdsim::samples::gpu

template class module<mdsim::samples::gpu::trajectory<3, double> >;
template class module<mdsim::samples::gpu::trajectory<2, double> >;

} // namespace halmd
