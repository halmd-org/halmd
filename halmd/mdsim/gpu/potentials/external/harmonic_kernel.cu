/*
 * Copyright © 2014 Felix Höfling
 * Copyright © 2020 Jaslo Ziska
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/mdsim/gpu/forces/external_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/external/harmonic_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {
namespace harmonic_kernel {

template <int dimension>
__device__ void harmonic<dimension>::fetch_param(unsigned int species)
{
    tie(offset_, stiffness_) <<= tex1Dfetch<float4>(t_param_, species);
}

} // namespace harmonic_kernel
} // namespace external
} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using namespace potentials::external::harmonic_kernel;

template class external_wrapper<3, harmonic<3> >;
template class external_wrapper<2, harmonic<2> >;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
