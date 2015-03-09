/*
 * Copyright © 2015 Nicolas Höft
 * Copyright © 2015 Felix Höfling
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

#define BOOST_TEST_MODULE copy_if

#include <algorithm>
#include <iomanip>
#include <random>
#include <stdint.h>

#include <boost/test/unit_test.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <halmd/algorithm/gpu/copy_if.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/cuda.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>
#include <test/unit/algorithm/gpu/copy_if_kernel.hpp>

BOOST_GLOBAL_FIXTURE( set_cuda_device )

/**
 * define test case template to test with various value types and predicates
 */
#define TEST_CASE_TYPE_PREDICATE(t)                         \
template <typename T, typename predicate_type>              \
void test_ ## t();                                          \
BOOST_AUTO_TEST_CASE(t)                                     \
{                                                           \
    test_ ## t<int8_t, select_all<int8_t>>();               \
    test_ ## t<int8_t, select_none<int8_t>>();              \
    test_ ## t<int8_t, select_even<int8_t>>();              \
    test_ ## t<int8_t, select_odd<int8_t>>();               \
    test_ ## t<int8_t, select_prime<int8_t>>();             \
    BOOST_TEST_MESSAGE("");                                 \
                                                            \
    test_ ## t<int, select_all<int>>();                     \
    test_ ## t<int, select_none<int>>();                    \
    test_ ## t<int, select_even<int>>();                    \
    test_ ## t<int, select_odd<int>>();                     \
    test_ ## t<int, select_prime<int>>();                   \
    BOOST_TEST_MESSAGE("");                                 \
                                                            \
    test_ ## t<int64_t, select_all<int64_t>>();             \
    test_ ## t<int64_t, select_none<int64_t>>();            \
    test_ ## t<int64_t, select_even<int64_t>>();            \
    test_ ## t<int64_t, select_odd<int64_t>>();             \
    test_ ## t<int64_t, select_prime<int64_t>>();           \
}                                                           \
template <typename T, typename predicate_type>              \
void test_ ## t()

TEST_CASE_TYPE_PREDICATE( copy_integer )
{
    cuda::host::vector<T> h_v(
        boost::make_counting_iterator(0)
      , boost::make_counting_iterator(12345679)
    );
    BOOST_TEST_MESSAGE("Processing array with " << h_v.size() << " elements using " <<
        halmd::demangled_name<predicate_type>()
    );

    cuda::vector<T> g_v_in(h_v.size());
    cuda::vector<T> g_v_out(h_v.size());
    cuda::copy(h_v.begin(), h_v.end(), g_v_in.begin());

    predicate_type pred;
    auto last_output = halmd::copy_if(g_v_in.begin(), g_v_in.end(), g_v_out.begin(), pred);
    std::size_t out_len = last_output - g_v_out.begin();

    BOOST_TEST_MESSAGE("CUDA copied " << out_len << " elements");
    cuda::host::vector<T> h_out(out_len);
    cuda::copy(g_v_out.begin(), g_v_out.begin() + out_len, h_out.begin());

    std::vector<T> h_reference(h_v.size());
    auto last_ref = std::copy_if(h_v.begin(), h_v.end(), h_reference.begin(), pred);
    h_reference.resize(last_ref - h_reference.begin());

    BOOST_TEST_MESSAGE( "CPU copied " << h_reference.size() << " elements");

    // compare results
    BOOST_CHECK_EQUAL_COLLECTIONS(
        h_out.begin()
      , h_out.end()
      , h_reference.begin()
      , h_reference.end()
    );
}

TEST_CASE_TYPE_PREDICATE( performance )
{
    BOOST_TEST_MESSAGE( "running " << halmd::demangled_name<predicate_type>() );

    cuda::host::vector<T> h_v(
        boost::make_counting_iterator(0)
      , boost::make_counting_iterator(12345679)
    );
    auto engine = std::default_random_engine{};
    std::shuffle(std::begin(h_v), std::end(h_v), engine);

    cuda::vector<T> g_v_in(h_v.size());
    cuda::vector<T> g_v_out(h_v.size());
    cuda::copy(h_v.begin(), h_v.end(), g_v_in.begin());

    predicate_type predicate;
    int const repeats = 5;
    double gpu_elapsed(0), host_elapsed(0);
    for (int i = 0; i < repeats; ++i) {
        auto fn = [&](double elapsed){gpu_elapsed += elapsed;};
        halmd::scoped_timer<halmd::timer> timer(fn);
        halmd::copy_if(g_v_in.begin(), g_v_in.end(), g_v_out.begin(), predicate);
    }
    std::vector<T> h_v_out(h_v.size());
    for (int i = 0; i < repeats; ++i) {
        auto fn = [&](double elapsed){host_elapsed += elapsed;};
        halmd::scoped_timer<halmd::timer> timer(fn);
        std::copy_if(h_v.begin(), h_v.end(), h_v_out.begin(), predicate);
    }
    host_elapsed /= repeats;
    gpu_elapsed /= repeats;
    double factor = host_elapsed / gpu_elapsed;
    BOOST_TEST_MESSAGE( std::setprecision(2) << std::fixed <<
        "host runtime: " << host_elapsed * 1e3 << " ms " <<
        "| gpu runtime: " << gpu_elapsed * 1e3 << " ms " <<
        "| factor: " << factor
    );

    if (host_elapsed > 1e-3) {
        BOOST_CHECK( factor >= 1 ); // select_none is optimised out by the host compiler
    }
}
