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

#include <cuda_wrapper.hpp>
#include "rand48.cuh"

namespace mdsim { namespace gpu { namespace rand48
{

extern cuda::function<void (uint48*)> leapfrog;
extern cuda::function<void (ushort3*, uint48 const*, uint48 const*, uint48*, uint48*, unsigned int)> set;
extern cuda::function<void (ushort3*, ushort3*)> save;
extern cuda::function<void (ushort3*, uint48 const*, uint48 const*, uint48*, uint48*, ushort3)> restore;
extern cuda::function<void (ushort3*, float*, unsigned int)> uniform;
extern cuda::function<void (ushort3*, uint*, unsigned int)> get;

extern cuda::symbol<uint48> a;
extern cuda::symbol<uint48> c;

template <typename T>
struct __boltzmann
{
    static cuda::function<void (T*, float, ushort3*)> ref;
};

template <typename T0, typename T1, typename T2>
inline void boltzmann(T0& t0, T1 const& t1, T2& t2)
{
    __boltzmann<typename T0::value_type>::ref(t0, t1, t2);
}

}}} // namespace mdsim::gpu::rand48

#endif /* ! MDSIM_GPU_RAND48_GLUE_HPP */
