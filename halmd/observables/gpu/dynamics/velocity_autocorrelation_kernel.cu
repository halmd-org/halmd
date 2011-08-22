/*
 * Copyright © 2008-2011  Peter Colberg
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

#include <halmd/observables/dynamics/velocity_autocorrelation.hpp>
#include <halmd/observables/gpu/dynamics/velocity_autocorrelation_kernel.hpp>
#include <halmd/observables/gpu/dynamics/tagged_particle.cuh>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

template <int dimension, unsigned int threads>
velocity_autocorrelation_wrapper<dimension, threads> const
velocity_autocorrelation_wrapper<dimension, threads>::wrapper = {
    accumulate<observables::dynamics::velocity_autocorrelation<vector_type>, threads>
};

// explicit template instantiation
template class velocity_autocorrelation_wrapper<3, 512>;
template class velocity_autocorrelation_wrapper<2, 512>;
template class velocity_autocorrelation_wrapper<3, 256>;
template class velocity_autocorrelation_wrapper<2, 256>;
template class velocity_autocorrelation_wrapper<3, 128>;
template class velocity_autocorrelation_wrapper<2, 128>;
template class velocity_autocorrelation_wrapper<3, 64>;
template class velocity_autocorrelation_wrapper<2, 64>;
template class velocity_autocorrelation_wrapper<3, 32>;
template class velocity_autocorrelation_wrapper<2, 32>;

} // namespace dynamics
} // namespace gpu
} // namespace observables
} // namespace halmd
