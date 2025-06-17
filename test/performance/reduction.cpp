/*
 * Copyright © 2021 Jaslo Ziska
 * Copyright © 2023 Felix Höfling
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

#define BOOST_TEST_MODULE reduction
#include <boost/test/unit_test.hpp>

#include <numeric>
#include <random>
#include <type_traits>

#include <boost/iterator/counting_iterator.hpp>

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/timer.hpp>
#include <test/performance/reduction_kernel.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>
#include <test/tools/dsfloat.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

BOOST_GLOBAL_FIXTURE(set_cuda_device);

template <typename T, typename K>
void reduce(cuda::memory::host::vector<T> const& h_input, cuda::memory::host::vector<T>& h_output, K& kernel)
{
    cuda::config dim(1, NTHREADS);

    cuda::memory::device::vector<T> g_input(h_input.size());
    cuda::memory::device::vector<T> g_output(h_output.size());

    cuda::event t1, t2;

    // copy data to gpu
    BOOST_CHECK(cuda::copy(
        h_input.begin()
      , h_input.end()
      , g_input.begin()) == g_input.end()
    );

    // configure kernel
    kernel.configure(dim.grid, dim.block);

    // record first event
    t1.record();
    // launch kernel
    kernel(g_input, g_output);
    cuda::thread::synchronize();
    // record second event
    t2.record();
    t2.synchronize();

    // copy result back to cpu
    BOOST_CHECK(cuda::copy(
        g_output.begin()
      , g_output.end()
      , h_output.begin()) == h_output.end()
    );

    float time = t2 - t1;
    BOOST_TEST_MESSAGE("summation of " << h_input.size() << " * " << halmd::demangled_name<T>() << ": " << 1000 * time << " ms");
}

/*
 * helper functions to compare collections of integers, floats, and fixed vectors
 */
// collection on integral type
template <typename Iter1, typename Iter2>
typename std::enable_if<std::is_integral<typename Iter1::value_type>::value, void>::type
check_equal_collections(Iter1 it1, Iter1 end1, Iter2 it2, Iter2 end2)
{
    BOOST_CHECK_EQUAL_COLLECTIONS(it1, end1, it2, end2);
}

// collection of floating-point type
template <typename Iter1, typename Iter2>
typename std::enable_if<std::is_floating_point<typename Iter1::value_type>::value, void>::type
check_equal_collections(Iter1 it1, Iter1 end1, Iter2 it2, Iter2 end2)
{
    typedef typename Iter1::value_type T;
    const T tolerance = 10 * dsfloat_aware_numeric_limits<T>::epsilon();
    while (it1 != end1 && it2 != end2) {
        BOOST_CHECK_CLOSE_FRACTION(*it1++, *it2++, tolerance);
    }
}

// collection on fixed_vector on integral type
template <typename Iter1, typename Iter2>
typename std::enable_if<std::is_enum<decltype(Iter1::value_type::static_size)>::value &&
                        std::is_integral<typename Iter1::value_type::value_type>::value, void>::type
check_equal_collections(Iter1 it1, Iter1 end1, Iter2 it2, Iter2 end2)
{
    BOOST_CHECK_EQUAL_COLLECTIONS(it1, end1, it2, end2);
}

// collection of fixed_vector on floating-point type
template <typename Iter1, typename Iter2>
typename std::enable_if<std::is_enum<decltype(Iter1::value_type::static_size)>::value &&
                        std::is_floating_point<typename Iter1::value_type::value_type>::value, void>::type
check_equal_collections(Iter1 it1, Iter1 end1, Iter2 it2, Iter2 end2)
{
    typedef typename Iter1::value_type::value_type T;
    static_assert(Iter1::value_type::static_size == Iter2::value_type::static_size, "mismatching sizes of array types");

    const T tolerance = 10 * dsfloat_aware_numeric_limits<T>::epsilon();
    while (it1 != end1 && it2 != end2) {
        // the following does not work due to fixed_vector<T>::begin returning T* instead of an iterator
        // check_equal_collections(it->begin(), it1->end(), it2->begin(), it2->end());
        for (unsigned int i = 0; i < it1->size(); ++i) {
            BOOST_CHECK_CLOSE_FRACTION((*it1)[i], (*it2)[i], tolerance);
        }
        ++it1, ++it2;
    }
}

/*
 * define test cases
 */
template <typename T, typename random_type = T>
void sum_integers_test()
{
    cuda::memory::host::vector<T> h_input(NTHREADS * NREDUCES);
    cuda::memory::host::vector<T> h_output(NREDUCES);
    std::vector<T> result(NREDUCES);

    std::default_random_engine gen;
    std::uniform_int_distribution<random_type> rand(0, 100);

    // generate random numbers
    std::generate(h_input.begin(), h_input.end(), std::bind(rand, std::ref(gen)));

    reduce(h_input, h_output, reduce_kernel<T>::kernel.sum);

    // calculate the result on the cpu
    for (unsigned int i = 0; i < NREDUCES; ++i) {
        result[i] = std::accumulate(h_input.begin() + NTHREADS * i, h_input.begin() + NTHREADS * (i + 1), T(0));
    }

    check_equal_collections(result.begin(), result.end(), h_output.begin(), h_output.end());
}

template <typename T, typename random_type = T>
void sum_reals_test()
{
    cuda::memory::host::vector<T> h_input(NTHREADS * NREDUCES);
    cuda::memory::host::vector<T> h_output(NREDUCES);
    std::vector<T> result(NREDUCES);

    std::default_random_engine gen;
    std::uniform_real_distribution<random_type> rand(0, 1);

    // generate random numbers
    std::generate(h_input.begin(), h_input.end(), std::bind(rand, std::ref(gen)));

    reduce(h_input, h_output, reduce_kernel<T>::kernel.sum);

    // calculate the result on the cpu
    for (unsigned int i = 0; i < NREDUCES; ++i) {
        result[i] = std::accumulate(h_input.begin() + NTHREADS * i, h_input.begin() + NTHREADS * (i + 1), T(0));
    }

    check_equal_collections(result.begin(), result.end(), h_output.begin(), h_output.end());
}

BOOST_AUTO_TEST_CASE(type_int)
{
    sum_integers_test<int>();
}

BOOST_AUTO_TEST_CASE(type_fixed_vector_int_2)
{
    sum_integers_test<halmd::fixed_vector<int, 2>, int>();
}

BOOST_AUTO_TEST_CASE(type_fixed_vector_int_3)
{
    sum_integers_test<halmd::fixed_vector<int, 3>, int>();
}

BOOST_AUTO_TEST_CASE(type_float)
{
    sum_reals_test<float>();
}

BOOST_AUTO_TEST_CASE(type_fixed_vector_float_2)
{
    sum_reals_test<halmd::fixed_vector<float, 2>, float>();
}

BOOST_AUTO_TEST_CASE(type_fixed_vector_float_3)
{
    sum_reals_test<halmd::fixed_vector<float, 3>, float>();
}

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
BOOST_AUTO_TEST_CASE(type_dsfloat)
{
    sum_reals_test<halmd::dsfloat, double>();
}

BOOST_AUTO_TEST_CASE(type_fixed_vector_dsfloat_3)
{
    sum_reals_test<halmd::fixed_vector<halmd::dsfloat, 3>, double>();
}
#endif

#ifdef USE_GPU_DOUBLE_PRECISION
BOOST_AUTO_TEST_CASE(type_double)
{
    sum_reals_test<double>();
}

BOOST_AUTO_TEST_CASE(type_fixed_vector_double_3)
{
    sum_reals_test<halmd::fixed_vector<double, 3>, double>();
}
#endif

