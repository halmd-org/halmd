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

#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/sample/trajectory.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
#endif
#include <halmd/mdsim/host/sample/trajectory.hpp>
#include <halmd/mdsim/samples/host/trajectory.hpp>
#include <halmd/util/logger.hpp>

namespace halmd
{
namespace mdsim { namespace samples { namespace host
{

template <int dimension>
typename trajectory<dimension, double>::pointer
trajectory<dimension, double>::create(options const& vm)
{
    return pointer(new mdsim::host::sample::trajectory<dimension, double>(vm));
}

template <int dimension>
typename trajectory<dimension, float>::pointer
trajectory<dimension, float>::create(options const& vm)
{
#ifdef WITH_CUDA
    if (module<mdsim::gpu::particle<dimension, float> >::fetch(vm)) {
        return pointer(new mdsim::gpu::sample::trajectory<trajectory<dimension, float> >(vm));
    }
#endif
    return pointer(new mdsim::host::sample::trajectory<dimension, float>(vm));
}

template class trajectory<3, double>;
template class trajectory<2, double>;
template class trajectory<3, float>;
template class trajectory<2, float>;

}}} // namespace mdsim::samples::host

template class module<mdsim::samples::host::trajectory<3, double> >;
template class module<mdsim::samples::host::trajectory<2, double> >;
template class module<mdsim::samples::host::trajectory<3, float> >;
template class module<mdsim::samples::host::trajectory<2, float> >;

} // namespace halmd
