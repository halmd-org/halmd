/*
 * Copyright © 2016 Daniel Kirchner
 * Copyright © 2017 Felix Höfling
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

#include <test/unit/numeric/mp/dsfloat_cuda_vector_kernel.hpp>

using namespace halmd;

namespace kernel {

__global__ void test_float4_ptr(float4 *g, fixed_vector<dsfloat, 3> increment)
{
    fixed_vector<dsfloat, 3> position;
    unsigned int species;

    tie (position, species) <<= tie(g[GTID], g[GTID + GTDIM]);

    position += increment;
    species *= 2;

    tie(g[GTID], g[GTID + GTDIM]) <<= tie(position, species);
}

__global__ void test_dsfloat_ptr(dsfloat_ptr<float4> g, fixed_vector<dsfloat, 3> increment)
{
    fixed_vector<dsfloat, 3> position;
    unsigned int species;

    tie (position, species) <<= g[GTID];

    position += increment;
    species *= 2;

    g[GTID] <<= tie(position, species);
}

template <typename T>
__global__ void overloaded_test(typename dsfloat_traits<T, 3>::ptr_type g, fixed_vector<T, 3> increment)
{
    fixed_vector<T, 3> value;
    int ignored;
    tie(value, ignored) <<= g[GTID];
    value += increment;
    g[GTID] <<= tie(value, ignored);
}

} // namespace kernel

dsfloat_kernel_wrapper dsfloat_kernel_wrapper::kernel = { kernel::test_float4_ptr, kernel::test_dsfloat_ptr };

template <typename float_type>
dsfloat_kernel_overloaded_wrapper<float_type> dsfloat_kernel_overloaded_wrapper<float_type>::kernel = {
    kernel::overloaded_test<float_type>
};

template class dsfloat_kernel_overloaded_wrapper<float>;
template class dsfloat_kernel_overloaded_wrapper<dsfloat>;
