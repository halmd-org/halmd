/*
 * Copyright © 2014 Felix Höfling
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

#define BOOST_TEST_MODULE modulation
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/utility/enable_if.hpp>
#include <limits>
#include <vector>

#include <halmd/observables/modulation.hpp>
#ifdef WITH_CUDA
# include <halmd/utility/gpu/device.hpp>
# include <test/unit/observables/modulation_kernel.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace boost::assign; // list_of, tuple_list_of
using namespace halmd;
using namespace halmd::observables::modulation;

/**
 * test modulation functions for density modes
 */

template <bool gpu, typename modulation_type, typename results_type>
typename boost::disable_if_c<gpu>::type
test_modulation(modulation_type const& modulation, results_type const& results)
{
    typedef typename modulation_type::argument_type vector_type;
    typedef typename modulation_type::result_type float_type;
    enum { dimension = vector_type::static_size };

    const float_type tolerance = std::max(float_type(1e-10), 2 * std::numeric_limits<float_type>::epsilon());

    BOOST_TEST_MESSAGE(modulation.message() << " on host");

    // evaluate modulation function at given set of points
    BOOST_FOREACH (typename results_type::value_type const& a, results) {
        vector_type r(0);
        r[dimension - 1] = get<0>(a);
        float_type f_r = modulation(r);
        BOOST_TEST_MESSAGE("r = (" << r << "), f(r) = " << f_r << ", expected: " << get<1>(a));
        BOOST_CHECK_CLOSE_FRACTION(f_r, get<1>(a), tolerance);
    }
}

#ifdef WITH_CUDA
template <bool gpu, typename modulation_type, typename results_type>
typename boost::enable_if_c<gpu>::type
test_modulation(modulation_type const& modulation, results_type const& results)
{
    typedef typename modulation_type::argument_type vector_type;
    typedef typename modulation_type::result_type float_type;
    enum { dimension = vector_type::static_size };

    typedef modulation_wrapper<dimension, modulation_type> wrapper_type;
    typedef typename wrapper_type::coalesced_vector_type coalesced_vector_type;

    const float_type tolerance = std::max(float_type(1e-10), 2 * std::numeric_limits<float_type>::epsilon());

    BOOST_TEST_MESSAGE(modulation.message() << " on GPU");

    // construct position vectors and copy to GPU
    unsigned int size = results.size();
    cuda::host::vector<coalesced_vector_type> h_r(size);
    cuda::host::vector<float_type> h_out(size);

    cuda::vector<coalesced_vector_type> g_r(h_r.size());
    cuda::vector<float_type> g_out(h_out.size());

    for (unsigned int i = 0; i < size; ++i) {
        vector_type r(0);
        r[dimension - 1] = get<0>(results[i]);
        h_r[i] = r;
    }
    cuda::copy(h_r, g_r);

    // call CUDA kernel
    cuda::configure(1, 32);     // 1 block of 32 threads
    wrapper_type::kernel.compute(g_r, g_out, size, modulation);
    cuda::thread::synchronize();

    // copy results to host and check against host implementation
    cuda::copy(g_out, h_out);
    for (unsigned int i = 0; i < size; ++i) {
        vector_type r = h_r[i];
        float_type f_r = modulation(r);
        BOOST_TEST_MESSAGE("r = (" << r << "), f(r) = " << h_out[i] << ", expected: " << f_r);
        BOOST_CHECK_CLOSE_FRACTION(h_out[i], f_r, tolerance);
    }
}
#endif

void test_unity()
{
    unity<1, double> modulation;

    BOOST_CHECK_EQUAL(modulation.message(), "");

    // evaluate some points of the modulation function
    std::vector<double> samples = list_of(-25.)(-1.)(0.)(0.1)(1.)(29);

    BOOST_FOREACH (double x, samples) {
        BOOST_CHECK_EQUAL(modulation(x), 1);
    }
}

template <int dimension, typename float_type, bool gpu>
void test_exponential()
{
    typedef exponential<dimension, float_type> modulation_type;

    // parameters: kappa, offset, box height
    float_type kappa = 0.1;
    float_type offset = 0;
    float_type box_height = 20;

    // expected results (z, f(r)) where r = (x, y, z) % L
    typedef boost::tuple<float_type, float_type> tuple_type;
    std::vector<tuple_type> results = tuple_list_of
        (-5. + box_height, 1.)
        (-1., 1.)
        (0., 1.)
        (0.1, 0.990049833749168)
        (1., 0.90483741803596)
        (9. - box_height, 0.406569659740599);    // test reduction to periodic box (-L/2, L/2)

    test_modulation<gpu>(modulation_type(kappa, offset, box_height), results);

    // second set of parameters, with negative sign of kappa and an offset
    kappa = -2;
    offset = 1;

    // expected results (z, f(r)) where r = (x, y, z) % L
    results = tuple_list_of
        (-5. - box_height, 6.14421235334284e-06)
        (-1., 0.0183156388887342)
        (0., 0.135335283236613)
        (0.1, 0.165298888221587)
        (1., 1.)
        (9. + box_height, 1.);

    test_modulation<gpu>(modulation_type(kappa, offset, box_height), results);
}

template <int dimension, typename float_type, bool gpu>
void test_exponential_slab()
{
    typedef exponential_slab<dimension, float_type> modulation_type;

    // parameters: kappa, width, offset, box height
    float_type kappa = 0.1;
    float_type width = 9;
    float_type offset = 0;
    float_type box_height = 20;

    // evaluate some points of the modulation function
    typedef boost::tuple<float_type, float_type> tuple_type;

    // expected results (z, f(r)) where r = (x, y, z) % L
    std::vector<tuple_type> results = tuple_list_of
        (-5. + box_height, 1.)
        (-1., 1.)
        (0., 1.)
        (0.1, 0.990049833749168)
        (1., 0.90483741803596)
        (9. - box_height, 0.);                   // test reduction to periodic box (-L/2, L/2)

    test_modulation<gpu>(modulation_type(kappa, width, offset, box_height), results);

    // second set of parameters, with negative sign of kappa and an offset
    kappa = -2;
    width = 2;
    offset = 1;

    // expected results (z, f(r)) where r = (x, y, z) % L
    results = tuple_list_of
        (-5. - box_height, 0.)
        (-1., 0.)
        (0., 0.135335283236613)
        (0.1, 0.165298888221587)
        (1., 1.)
        (9. + box_height, 1.);

    test_modulation<gpu>(modulation_type(kappa, width, offset, box_height), results);
}

template <int dimension, typename float_type, bool gpu>
void test_catenary()
{
    typedef catenary<dimension, float_type> modulation_type;

    // parameters: kappa, width, offset, box height
    float_type kappa = 0.1;
    float_type width = 10;
    float_type offset = 0;
    float_type box_height = 20;

    // evaluate some points of the modulation function
    typedef boost::tuple<float_type, float_type> tuple_type;

    // expected results (z, f(r)) where r = (x, y, z) % L
    std::vector<tuple_type> results = tuple_list_of
        (-5. + box_height, 1.)
        (-1., 0.89125667470052)
        (0., 0.886818883970074)
        (0.1, 0.886863225283782)
        (1., 0.89125667470052)
        (9. - box_height, 1.);                   // test reduction to periodic box (-L/2, L/2)

    test_modulation<gpu>(modulation_type(kappa, width, offset, box_height), results);

    // second set of parameters, with negative sign of kappa and an offset
    kappa = -2;
    width = 10;
    offset = 5;

    // expected results (z, f(r)) where r = (x, y, z) % L
    results = tuple_list_of
        (-5. - box_height, 1.)
        (-1., 1.)
        (0., 1.)
        (0.1, 0.818730753907951)
        (1., 0.135335298187646)
        (9. + box_height, 0.135335298187646);

    test_modulation<gpu>(modulation_type(kappa, width, offset, box_height), results);
}

BOOST_AUTO_TEST_CASE( modulation_host_unity ) {
    test_unity();
}

BOOST_AUTO_TEST_CASE( modulation_host_exponential ) {
    test_exponential<2, float, false>();
    test_exponential<3, float, false>();
    test_exponential<2, double, false>();
    test_exponential<3, double, false>();
}

BOOST_AUTO_TEST_CASE( modulation_host_exponential_slab ) {
    test_exponential_slab<2, float, false>();
    test_exponential_slab<3, float, false>();
    test_exponential_slab<2, double, false>();
    test_exponential_slab<3, double, false>();
}

BOOST_AUTO_TEST_CASE( modulation_host_catenary ) {
    test_catenary<2, float, false>();
    test_catenary<3, float, false>();
    test_catenary<2, double, false>();
    test_catenary<3, double, false>();
}

#ifdef WITH_CUDA
BOOST_FIXTURE_TEST_CASE( modulation_gpu_exponential, device ) {
    test_exponential<2, float, true>();
    test_exponential<3, float, true>();
}

BOOST_FIXTURE_TEST_CASE( modulation_gpu_exponential_slab, device ) {
    test_exponential_slab<2, float, true>();
    test_exponential_slab<3, float, true>();
}

BOOST_FIXTURE_TEST_CASE( modulation_gpu_catenary, device ) {
    test_catenary<2, float, true>();
    test_catenary<3, float, true>();
}
#endif // WITH_CUDA

