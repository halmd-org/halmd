/* DSFUN addition test
 *
 * Copyright © 2009  Peter Colberg
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

#ifndef TEST_DSFUN_ADD_HPP
#define TEST_DSFUN_ADD_HPP

#include <cuda_wrapper.hpp>
#include <halmd/numeric/gpu/blas/dsfloat.cuh>

using halmd::numeric::gpu::blas::dsfloat;

extern cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_add;

#endif /* ! TEST_DSFUN_ADD_HPP */
