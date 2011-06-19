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

#ifndef HALMD_MDSIM_GPU_NEIGHBOUR_WITHOUT_BINNING_KERNEL_HPP
#define HALMD_MDSIM_GPU_NEIGHBOUR_WITHOUT_BINNING_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension>
struct neighbour_without_binning_wrapper
{
    /** (cutoff lengths + neighbour list skin)² */
    cuda::texture<float> rr_cut_skin;
    /** update neighbour lists */
    cuda::function<void (
        float4 const*
      , unsigned int
      , float4 const*
      , unsigned int
      , fixed_vector<float, dimension>
      , unsigned int*
      , unsigned int
      , unsigned int
      , int*
    )> update;

    static neighbour_without_binning_wrapper kernel;
};

template <int dimension>
neighbour_without_binning_wrapper<dimension> const& get_neighbour_without_binning_kernel()
{
    return neighbour_without_binning_wrapper<dimension>::kernel;
}

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOUR_WITHOUT_BINNING_KERNEL_HPP */
