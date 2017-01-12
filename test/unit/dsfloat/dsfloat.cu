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

#include <test/unit/dsfloat/dsfloat.hpp>

using namespace halmd;

__global__ void dsfloat_kernel_test_1 (float4 *g, fixed_vector<dsfloat,3> increment) {
    fixed_vector<dsfloat,3> position;
    unsigned int species;

    tie (position, species) <<= tie(g[GTID], g[GTID + GTDIM]);

    for (int i = 0; i < 1000; i++) {
        position += increment;
    }

    tie(g[GTID], g[GTID + GTDIM]) <<= tie(position, species);
}

__global__ void dsfloat_kernel_test_2 (dsfloat_ptr<float4> g, fixed_vector<dsfloat,3> increment) {
    fixed_vector<dsfloat,3> position;
    unsigned int species;

    tie (position, species) <<= g[GTID];

    for (int i = 0; i < 1000; i++) {
        position += increment;
    }

    g[GTID] <<= tie(position, species);
}

dsfloat_kernel_wrapper dsfloat_kernel_wrapper::kernel = {
        dsfloat_kernel_test_1,
        dsfloat_kernel_test_2
};

template<typename T>
__global__ void dsfloat_kernel_overloaded_test (typename dsfloat_traits<T, 3>::ptr_type g, fixed_vector<T,3> increment)
{
    fixed_vector<T, 3> value;
    int ignored;
    tie(value, ignored) <<= g[GTID];
    value += increment;
    g[GTID] <<= tie(value, ignored);
}

template<typename float_type>
dsfloat_kernel_overloaded_wrapper<float_type> dsfloat_kernel_overloaded_wrapper<float_type>::kernel = {
        dsfloat_kernel_overloaded_test<float_type>
};

template class dsfloat_kernel_overloaded_wrapper<float>;
template class dsfloat_kernel_overloaded_wrapper<dsfloat>;
