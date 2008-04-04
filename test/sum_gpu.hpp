/* Test alternatives for summing over a CUDA device vector
 *
 * Copyright (C) 2008  Peter Colberg
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

#ifndef MDSIM_TEST_SUMAVG
#define MDSIM_TEST_SUMAVG

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace mdsim { namespace test { namespace gpu {

extern cuda::function<void (float*, float*)> blocksum;

}}} // namespace mdsim::test::gpu

#endif /* ! MDSIM_TEST_SUMAVG */
