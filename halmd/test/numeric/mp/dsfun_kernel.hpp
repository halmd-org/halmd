/*
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

#ifndef HALMD_TEST_NUMERIC_MP_DSFUN_HPP
#define HALMD_TEST_NUMERIC_MP_DSFUN_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

namespace halmd
{
namespace test
{

extern cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_add;
extern cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_sub;
extern cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_mul;
extern cuda::function<void (float const*, float const*, dsfloat*)> kernel_mulss;
extern cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_div;
extern cuda::function<void (dsfloat const*, dsfloat*)> kernel_sqrt;

} // namespace test

} // namespace halmd

#endif /* ! HALMD_TEST_NUMERIC_MP_DSFUN_HPP */
