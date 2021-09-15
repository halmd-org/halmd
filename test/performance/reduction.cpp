/*
 * Copyright © 2021 Jaslo Ziska
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

#include <boost/iterator/counting_iterator.hpp>

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/timer.hpp>
#include <test/performance/reduction_kernel.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>

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
    BOOST_TEST_MESSAGE("summation of " << h_input.size() << " * " << halmd::demangled_name<T>() << ": " << time << " s");
}

BOOST_AUTO_TEST_CASE(type_int)
{
    cuda::memory::host::vector<int> h_input(NTHREADS * NREDUCES);
    cuda::memory::host::vector<int> h_output(NREDUCES);
    std::vector<int> result(NREDUCES);

    std::default_random_engine gen;
    std::uniform_int_distribution<int> rand(0, 100);

    // generate random numbers
    std::generate(h_input.begin(), h_input.end(), std::bind(rand, std::ref(gen)));

    reduce(h_input, h_output, reduce_int_kernel);

    // calculate the result on the cpu
    for (unsigned int i = 0; i < NREDUCES; ++i) {
        result[i] = std::accumulate(h_input.begin() + NTHREADS * i, h_input.begin() + NTHREADS * (i + 1), 0.0);
    }

    BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(), h_output.begin(), h_output.end());
}

#ifdef USE_GPU_SINGLE_PRECISION
BOOST_AUTO_TEST_CASE(type_float)
{
    cuda::memory::host::vector<float> h_input(NTHREADS * NREDUCES);
    cuda::memory::host::vector<float> h_output(NREDUCES);
    std::vector<float> result(NREDUCES);

    std::default_random_engine gen;
    std::uniform_real_distribution<float> rand(0, 1);

    // generate random numbers
    std::generate(h_input.begin(), h_input.end(), std::bind(rand, std::ref(gen)));

    reduce(h_input, h_output, reduce_float_kernel);

    // calculate the result on the cpu
    for (unsigned int i = 0; i < NREDUCES; ++i) {
        result[i] = std::accumulate(h_input.begin() + NTHREADS * i, h_input.begin() + NTHREADS * (i + 1), 0.0);
    }

    for (unsigned int i = 0; i < NREDUCES; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(result[i], h_output[i], 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(type_fixed_vector_float)
{
    using fixed_vector_type = halmd::fixed_vector<float, 3>;

    cuda::memory::host::vector<fixed_vector_type> h_input(NTHREADS * NREDUCES);
    cuda::memory::host::vector<fixed_vector_type> h_output(NREDUCES);
    std::vector<fixed_vector_type> result(NREDUCES);

    std::default_random_engine gen;
    std::uniform_real_distribution<float> rand(0, 1);

    // generate random numbers
    std::generate(h_input.begin(), h_input.end(), std::bind(rand, std::ref(gen)));

    reduce(h_input, h_output, reduce_fixed_vector_float_kernel);

    // calculate the result on the cpu
    for (unsigned int i = 0; i < NREDUCES; ++i) {
        result[i] = std::accumulate(
            h_input.begin() + NTHREADS * i
          , h_input.begin() + NTHREADS * (i + 1)
          , fixed_vector_type(0.0)
        );
    }

    for (unsigned int i = 0; i < NREDUCES; ++i) {
        for (unsigned int j = 0; j < fixed_vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION(result[i][j], h_output[i][j], 1e-6);
        }
    }
}
#endif

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
BOOST_AUTO_TEST_CASE(type_dsfloat)
{
    using halmd::dsfloat;

    cuda::memory::host::vector<dsfloat> h_input(NTHREADS * NREDUCES);
    cuda::memory::host::vector<dsfloat> h_output(NREDUCES);
    std::vector<double> result(NREDUCES);

    std::default_random_engine gen;
    std::uniform_real_distribution<double> rand(0, 1);

    // generate random numbers
    std::generate(h_input.begin(), h_input.end(), std::bind(rand, std::ref(gen)));

    reduce(h_input, h_output, reduce_dsfloat_kernel);

    // calculate the result on the cpu
    for (unsigned int i = 0; i < NREDUCES; ++i) {
        result[i] = std::accumulate(h_input.begin() + NTHREADS * i, h_input.begin() + NTHREADS * (i + 1), 0.0);
    }

    std::vector<double> h_result(h_output.begin(), h_output.end());
    for (unsigned int i = 0; i < NREDUCES; ++i) {
        BOOST_CHECK_CLOSE_FRACTION(result[i], h_result[i], 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(type_fixed_vector_dsfloat)
{
    using halmd::dsfloat;
    using fixed_vector_type = halmd::fixed_vector<dsfloat, 3>;

    cuda::memory::host::vector<fixed_vector_type> h_input(NTHREADS * NREDUCES);
    cuda::memory::host::vector<fixed_vector_type> h_output(NREDUCES);
    std::vector<fixed_vector_type> result(NREDUCES);

    std::default_random_engine gen;
    std::uniform_real_distribution<double> rand(0, 1);

    // generate random numbers
    std::generate(h_input.begin(), h_input.end(), std::bind(rand, std::ref(gen)));

    reduce(h_input, h_output, reduce_fixed_vector_dsfloat_kernel);

    // calculate the result on the cpu
    for (unsigned int i = 0; i < NREDUCES; ++i) {
        result[i] = std::accumulate(
            h_input.begin() + NTHREADS * i
          , h_input.begin() + NTHREADS * (i + 1)
          , fixed_vector_type(0.0)
        );
    }

    for (unsigned int i = 0; i < NREDUCES; ++i) {
        for (unsigned int j = 0; j < fixed_vector_type::static_size; ++j) {
            BOOST_CHECK_CLOSE_FRACTION(result[i][j], h_output[i][j], 1e-12);
        }
    }
}
#endif
