/*
 * Copyright © 2016 Daniel Kirchner
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

#ifndef HALMD_TEST_UNIT_DSFLOAT_DSFLOAT_HPP
#define HALMD_TEST_UNIT_DSFLOAT_DSFLOAT_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

struct dsfloat_kernel_wrapper
{
    cuda::function <void (float4*, halmd::fixed_vector<halmd::dsfloat,3>)> test1;
    cuda::function <void (halmd::dsfloat_ptr<float4>, halmd::fixed_vector<halmd::dsfloat,3>)> test2;

    static dsfloat_kernel_wrapper kernel;
};



template<typename value_type, int dimension>
struct dsfloat_traits {
    typedef typename halmd::mdsim::type_traits<dimension, value_type>::gpu::coalesced_vector_type* ptr_type;
};

template<int dimension>
struct dsfloat_traits<halmd::dsfloat, dimension> {
    typedef halmd::dsfloat_ptr<typename halmd::mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type> ptr_type;
};

template<typename float_type>
struct dsfloat_kernel_overloaded_wrapper
{
    typedef typename dsfloat_traits<float_type, 3>::ptr_type ptr_type;
    cuda::function <void (ptr_type, halmd::fixed_vector<float_type,3>)> test;

    static dsfloat_kernel_overloaded_wrapper kernel;
};

#endif /* ! HALMD_TEST_UNIT_DSFLOAT_DSFLOAT_HPP */
