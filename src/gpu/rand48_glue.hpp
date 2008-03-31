/* Parallelized rand48 random number generator for CUDA
 *
 * Copyright (C) 2007  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef MDSIM_GPU_RAND48_GLUE_HPP
#define MDSIM_GPU_RAND48_GLUE_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace mdsim { namespace gpu
{

template <typename T>
class rand48
{
public:
    static cuda::symbol<uint3> a;
    static cuda::symbol<uint3> c;

    static cuda::function<void (ushort3*, uint3*, uint3*, unsigned int)> init;
    static cuda::function<void (ushort3*, ushort3*)> save;
    static cuda::function<void (ushort3*, uint3*, uint3*, ushort3)> restore;
    static cuda::function<void (ushort3*, float*, unsigned int)> get_uniform;
    static cuda::function<void (ushort3*, T*, unsigned int)> get_unit_vector;
};

}} // namespace mdsim::gpu

#endif /* ! MDSIM_GPU_RAND48_GLUE_HPP */
