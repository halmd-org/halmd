/*
 * Copyright © 2012 Peter Colberg
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

#include <halmd/algorithm/gpu/reduce_kernel.cuh>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <test/unit/algorithm/gpu/reduce_kernel.hpp>

template class halmd::reduction_kernel<sum<float, halmd::dsfloat> >;
template class halmd::reduction_kernel<sum_with_constant<float, halmd::dsfloat> >;
template class halmd::reduction_kernel<sum_of_squares<float, halmd::dsfloat> >;
template class halmd::reduction_kernel<sum_of_squares_with_constant<float, halmd::dsfloat> >;
template class halmd::reduction_kernel<sum_of_cubes<float, halmd::dsfloat> >;
#ifdef HALMD_GPU_DOUBLE_PRECISION
template class halmd::reduction_kernel<sum<float, double> >;
template class halmd::reduction_kernel<sum_with_constant<float, double> >;
template class halmd::reduction_kernel<sum_of_squares<float, double> >;
template class halmd::reduction_kernel<sum_of_squares_with_constant<float, double> >;
template class halmd::reduction_kernel<sum_of_cubes<float, double> >;
#endif
