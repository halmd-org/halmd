/*
 * Copyright © 2008, 2012 Peter Colberg
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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE radix_sort
#include <boost/test/unit_test.hpp>

#include <halmd/algorithm/host/radix_sort.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/algorithm/gpu/radix_sort.hpp>
# include <test/tools/cuda.hpp>
#endif

#include <boost/iterator/transform_iterator.hpp>

#include <algorithm>
#include <random>
#include <vector>

/**
 * Returns array of uniformly distributed unsigned integers.
 */
static std::vector<unsigned int> make_uniform_array(int count)
{
    std::vector<unsigned int> output(count);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> dist;
    std::generate(output.begin(), output.end(), [&]() {
        return dist(gen);
    });
    return std::move(output);
}

/**
 * Test std::sort.
 */
static void test_std_sort(int count, int repeat)
{
    std::vector<unsigned int> input = make_uniform_array(count);

    BOOST_TEST_MESSAGE( "  " << count << " elements" );
    BOOST_TEST_MESSAGE( "  " << repeat << " iterations" );

    halmd::accumulator<double> elapsed;
    for (int i = 0; i < repeat; ++i) {
        std::vector<unsigned int> output(input.begin(), input.end());
        {
            halmd::scoped_timer<halmd::timer> t(elapsed);
            std::sort(output.begin(), output.end());
        }
        BOOST_CHECK( std::is_sorted(output.begin(), output.end()) );
    }
    BOOST_TEST_MESSAGE( "  " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration" );
}

/**
 * Test halmd::radix_sort on host.
 */
static void test_radix_sort_host(int count, int repeat)
{
    std::vector<unsigned int> input = make_uniform_array(count);
    std::vector<unsigned int> result(input.begin(), input.end());
    std::sort(result.begin(), result.end());

    BOOST_TEST_MESSAGE( "  " << count << " elements" );
    BOOST_TEST_MESSAGE( "  " << repeat << " iterations" );

    halmd::accumulator<double> elapsed;
    for (int i = 0; i < repeat; ++i) {
        std::vector<unsigned int> output(input.begin(), input.end());
        {
            halmd::scoped_timer<halmd::timer> t(elapsed);
            halmd::radix_sort(output.begin(), output.end());
        }
        BOOST_CHECK_EQUAL_COLLECTIONS(
            output.begin()
          , output.end()
          , result.begin()
          , result.end()
        );
    }
    BOOST_TEST_MESSAGE( "  " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration" );
}

#ifdef HALMD_WITH_GPU
/**
 * Test halmd::radix_sort on GPU.
 */
static void test_radix_sort_gpu(int count, int repeat, int threads)
{
    std::vector<unsigned int> input = make_uniform_array(count);
    cuda::vector<unsigned int> g_input(count);
    BOOST_CHECK( cuda::copy(
        input.begin()
      , input.end()
      , g_input.begin()) == g_input.end()
    );
    std::vector<unsigned int> result(input.begin(), input.end());
    std::sort(result.begin(), result.end());

    BOOST_TEST_MESSAGE( "  " << count << " elements" );
    BOOST_TEST_MESSAGE( "  " << repeat << " iterations" );
    BOOST_TEST_MESSAGE( "  " << threads << " threads per block" );

    halmd::accumulator<double> elapsed;
    halmd::algorithm::gpu::radix_sort<unsigned int> radix_sort(count, threads);
    for (int i = 0; i < repeat; ++i) {
        cuda::vector<unsigned int> g_output(count);
        cuda::vector<unsigned int> g_value(count);
        BOOST_CHECK( cuda::copy(
            g_input.begin()
          , g_input.end()
          , g_output.begin()) == g_output.end()
        );
        {
            halmd::scoped_timer<halmd::timer> t(elapsed);
            radix_sort(g_output, g_value);
            cuda::thread::synchronize();
        }
        cuda::host::vector<unsigned int> h_output(count);
        BOOST_CHECK( cuda::copy(
            g_output.begin()
          , g_output.end()
          , h_output.begin()) == h_output.end()
        );
        BOOST_CHECK_EQUAL_COLLECTIONS(
            h_output.begin()
          , h_output.end()
          , result.begin()
          , result.end()
        );
    }
    BOOST_TEST_MESSAGE( "  " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration" );
}

/**
 * Test generation of permutation using halmd::radix_sort on GPU.
 */
static void test_permutation_gpu(int count, int repeat, int threads)
{
    std::vector<unsigned int> input_key = make_uniform_array(count);
    cuda::vector<unsigned int> g_input_key(count);
    BOOST_CHECK( cuda::copy(
        input_key.begin()
      , input_key.end()
      , g_input_key.begin()) == g_input_key.end()
    );
    std::vector<unsigned int> result(input_key.begin(), input_key.end());
    std::sort(result.begin(), result.end());

    std::vector<unsigned int> input_value(count);
    std::iota(input_value.begin(), input_value.end(), 0);
    cuda::vector<unsigned int> g_input_value(count);
    BOOST_CHECK( cuda::copy(
        input_value.begin()
      , input_value.end()
      , g_input_value.begin()) == g_input_value.end()
    );

    BOOST_TEST_MESSAGE( "  " << count << " elements" );
    BOOST_TEST_MESSAGE( "  " << repeat << " iterations" );
    BOOST_TEST_MESSAGE( "  " << threads << " threads per block" );

    halmd::accumulator<double> elapsed;
    halmd::algorithm::gpu::radix_sort<unsigned int> radix_sort(count, threads);
    for (int i = 0; i < repeat; ++i) {
        cuda::vector<unsigned int> g_output_key(count);
        BOOST_CHECK( cuda::copy(
            g_input_key.begin()
          , g_input_key.end()
          , g_output_key.begin()) == g_output_key.end()
        );

        cuda::vector<unsigned int> g_output_value(count);
        BOOST_CHECK( cuda::copy(
            g_input_value.begin()
          , g_input_value.end()
          , g_output_value.begin()) == g_output_value.end()
        );

        {
            halmd::scoped_timer<halmd::timer> t(elapsed);
            radix_sort(g_output_key, g_output_value);
            cuda::thread::synchronize();
        }

        cuda::host::vector<unsigned int> h_output_key(count);
        BOOST_CHECK( cuda::copy(
            g_output_key.begin()
          , g_output_key.end()
          , h_output_key.begin()) == h_output_key.end()
        );
        BOOST_CHECK_EQUAL_COLLECTIONS(
            h_output_key.begin()
          , h_output_key.end()
          , result.begin()
          , result.end()
        );

        cuda::host::vector<unsigned int> h_output_value(count);
        BOOST_CHECK( cuda::copy(
            g_output_value.begin()
          , g_output_value.end()
          , h_output_value.begin()) == h_output_value.end()
        );
        auto value_to_key = [&](unsigned int value) {
            return input_key[value];
        };
        BOOST_CHECK_EQUAL_COLLECTIONS(
            boost::make_transform_iterator(h_output_value.begin(), value_to_key)
          , boost::make_transform_iterator(h_output_value.end(), value_to_key)
          , result.begin()
          , result.end()
        );
    }
    BOOST_TEST_MESSAGE( "  " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration" );
}
#endif /* HALMD_WITH_GPU */

HALMD_TEST_INIT( radix_sort )
{
    using namespace boost::unit_test;

    for (int count : {1000000, 100000, 10000, 1024, 1000, 512, 100, 10, 2, 1}) {
        test_suite* ts = BOOST_TEST_SUITE( std::to_string(count) );
        framework::master_test_suite().add(ts);
        {
            int const repeat = std::max(100000 / count, 5);
            auto std_sort = [=]() {
                test_std_sort(count, repeat);
            };
            ts->add(BOOST_TEST_CASE( std_sort ));
        }
        {
            int const repeat = std::max(100000 / count, 5);
            auto radix_sort_host = [=]() {
                test_radix_sort_host(count, repeat);
            };
            ts->add(BOOST_TEST_CASE( radix_sort_host ));
        }
#ifdef HALMD_WITH_GPU
        for (int threads : {128, 256}) {
            int const repeat = std::max(100 / count, 5);
            auto radix_sort_gpu = [=]() {
                set_cuda_device device;
                test_radix_sort_gpu(count, repeat, threads);
            };
            ts->add(BOOST_TEST_CASE( radix_sort_gpu ));

            auto permutation_gpu = [=]() {
                set_cuda_device device;
                test_permutation_gpu(count, repeat, threads);
            };
            ts->add(BOOST_TEST_CASE( permutation_gpu ));
        }
#endif
    }
}
