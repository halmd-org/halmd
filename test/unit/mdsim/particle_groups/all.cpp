/*
 * Copyright © 2012 Peter Colberg
 * Copyright © 2020 Jaslo Ziska
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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE all
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_groups/all.hpp>
# include <test/tools/cuda.hpp>
#endif
#define HALMD_TEST_NO_LOGGING
#include <test/tools/ctest.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <algorithm>
#include <numeric>

/**
 * Returns randomly shuffled particle ids.
 */
template <typename size_type>
static std::vector<size_type>
make_random_id(size_type nparticle)
{
    std::vector<size_type> id(nparticle);
    std::iota(id.begin(), id.end(), 0);
    std::random_shuffle(id.begin(), id.end());
    return id;
}

/**
 * Test particle_group::ordered()
 */
template <typename test_suite_type>
static void
test_ordered(
    unsigned int nparticle
  , unsigned int nspecies
  , unsigned int repeat
)
{
    typedef typename test_suite_type::particle_type particle_type;
    typedef typename test_suite_type::particle_group_type particle_group_type;
    typedef typename particle_group_type::array_type array_type;
    typedef typename particle_group_type::size_type size_type;
    typedef typename test_suite_type::all_type all_type;

    BOOST_TEST_MESSAGE("  " << "all of " << nparticle << " particles");
    BOOST_TEST_MESSAGE("  " << repeat << " iterations");

    std::shared_ptr<particle_type> particle = std::make_shared<particle_type>(nparticle, nspecies);
    std::shared_ptr<particle_group_type> group = std::make_shared<all_type>(particle);
    {
        BOOST_CHECK_EQUAL( *group->size(), nparticle );
    }

    halmd::cache<> size_cache = group->size();
    halmd::cache<> ordered_cache = group->ordered();

    halmd::accumulator<double> elapsed;
    for (size_type i = 0; i < repeat; ++i) {
        std::vector<size_type> reverse_id = make_random_id(nparticle);

        BOOST_CHECK( size_cache == group->size() );
        BOOST_CHECK( ordered_cache == group->ordered() );
        BOOST_CHECK( set_reverse_id(
            *particle
          , reverse_id.begin()) == reverse_id.end()
        );
        {
            halmd::scoped_timer<halmd::timer> t(elapsed);
            halmd::cache<array_type> const& ordered = group->ordered();
            BOOST_CHECK( ordered_cache != ordered );
            ordered_cache = ordered;
        }

        std::vector<size_type> ordered(nparticle);
        BOOST_CHECK( get_ordered(
            *group
          , ordered.begin()) == ordered.end()
        );
        BOOST_CHECK_EQUAL_COLLECTIONS(
            ordered.begin()
          , ordered.end()
          , reverse_id.begin()
          , reverse_id.end()
        );
    }
    BOOST_TEST_MESSAGE( "  " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration" );
}

/**
 * Test particle_group::unordered()
 */
template <typename test_suite_type>
static void
test_unordered(
    unsigned int nparticle
  , unsigned int nspecies
  , unsigned int repeat
)
{
    typedef typename test_suite_type::particle_type particle_type;
    typedef typename test_suite_type::particle_group_type particle_group_type;
    typedef typename particle_group_type::size_type size_type;
    typedef typename test_suite_type::all_type all_type;

    BOOST_TEST_MESSAGE("  " << "all of " << nparticle << " particles");
    BOOST_TEST_MESSAGE("  " << repeat << " iterations");

    std::shared_ptr<particle_type> particle = std::make_shared<particle_type>(nparticle, nspecies);
    std::shared_ptr<particle_group_type> group = std::make_shared<all_type>(particle);
    {
        BOOST_CHECK_EQUAL( *group->size(), nparticle );
    }

    halmd::cache<> size_cache = group->size();
    halmd::cache<> unordered_cache = group->unordered();

    halmd::accumulator<double> elapsed;
    for (size_type i = 0; i < repeat; ++i) {
        std::vector<size_type> reverse_id = make_random_id(nparticle);

        BOOST_CHECK( size_cache == group->size() );
        BOOST_CHECK( unordered_cache == group->unordered() );
        BOOST_CHECK( set_reverse_id(
            *particle
          , reverse_id.begin()) == reverse_id.end()
        );
        {
            halmd::scoped_timer<halmd::timer> t(elapsed);
            BOOST_CHECK( unordered_cache == group->unordered() );
        }

        std::vector<size_type> unordered(nparticle);
        BOOST_CHECK( get_unordered(
            *group
          , unordered.begin()) == unordered.end()
        );
        BOOST_CHECK_EQUAL_COLLECTIONS(
            unordered.begin()
          , unordered.end()
          , boost::make_counting_iterator(size_type(0))
          , boost::make_counting_iterator(nparticle)
        );
    }
    BOOST_TEST_MESSAGE( "  " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration" );
}

/**
 * Classes for host test suite.
 */
template <int dimension, typename float_type>
struct test_suite_host
{
    typedef halmd::mdsim::host::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::host::particle_group particle_group_type;
    typedef halmd::mdsim::host::particle_groups::all<particle_type> all_type;
};

#ifdef HALMD_WITH_GPU
/**
 * Classes for GPU test suite.
 */
template <int dimension, typename float_type>
struct test_suite_gpu
{
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::gpu::particle_group particle_group_type;
    typedef halmd::mdsim::gpu::particle_groups::all<particle_type> all_type;
};
#endif

/**
 * Data-driven test case registration.
 */
using namespace boost::unit_test;

unsigned int const DATA_ARRAY[] = {500000, 25000, 1000};
auto dataset = data::make(DATA_ARRAY);
unsigned int const nspecies = 1;
unsigned int constexpr repeat = 10;

BOOST_AUTO_TEST_SUITE( host )
    BOOST_AUTO_TEST_SUITE( two )
#ifdef USE_HOST_SINGLE_PRECISION
        typedef test_suite_host<2, float> test_suite_type;
#else
        typedef test_suite_host<2, double> test_suite_type;
#endif
        BOOST_DATA_TEST_CASE( ordered, dataset, nparticle ) {
            test_ordered<test_suite_type>(nparticle, nspecies, repeat);
        }
        BOOST_DATA_TEST_CASE( unordered, dataset, nparticle ) {
            test_unordered<test_suite_type>(nparticle, nspecies, repeat);
        }
    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE( three )
#ifdef USE_HOST_SINGLE_PRECISION
        typedef test_suite_host<3, float> test_suite_type;
#else
        typedef test_suite_host<3, double> test_suite_type;
#endif
        BOOST_DATA_TEST_CASE( ordered, dataset, nparticle ) {
            test_ordered<test_suite_type>(nparticle, nspecies, repeat);
        }
        BOOST_DATA_TEST_CASE( unordered, dataset, nparticle ) {
            test_unordered<test_suite_type>(nparticle, nspecies, repeat);
        }
    BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#ifdef HALMD_WITH_GPU
BOOST_FIXTURE_TEST_SUITE( gpu, set_cuda_device )
    BOOST_AUTO_TEST_SUITE( two )
# ifdef USE_GPU_SINGLE_PRECISION
        BOOST_AUTO_TEST_SUITE( type_float )
            typedef test_suite_gpu<2, float> test_suite_type;

            BOOST_DATA_TEST_CASE( ordered, dataset, nparticle ) {
                test_ordered<test_suite_type>(nparticle, nspecies, repeat);
            }
            BOOST_DATA_TEST_CASE( unordered, dataset, nparticle ) {
                test_unordered<test_suite_type>(nparticle, nspecies, repeat);
            }
        BOOST_AUTO_TEST_SUITE_END()
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
        BOOST_AUTO_TEST_SUITE( type_dsfloat )
            typedef test_suite_gpu<2, halmd::dsfloat> test_suite_type;

            BOOST_DATA_TEST_CASE( ordered, dataset, nparticle ) {
                test_ordered<test_suite_type>(nparticle, nspecies, repeat);
            }
            BOOST_DATA_TEST_CASE( unordered, dataset, nparticle ) {
                test_unordered<test_suite_type>(nparticle, nspecies, repeat);
            }
        BOOST_AUTO_TEST_SUITE_END()
# endif
    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE( three )
# ifdef USE_GPU_SINGLE_PRECISION
        BOOST_AUTO_TEST_SUITE( type_float )
            typedef test_suite_gpu<3, float> test_suite_type;

            BOOST_DATA_TEST_CASE( ordered, dataset, nparticle ) {
                test_ordered<test_suite_type>(nparticle, nspecies, repeat);
            }
            BOOST_DATA_TEST_CASE( unordered, dataset, nparticle ) {
                test_unordered<test_suite_type>(nparticle, nspecies, repeat);
            }
        BOOST_AUTO_TEST_SUITE_END()
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
        BOOST_AUTO_TEST_SUITE( type_dsfloat )
            typedef test_suite_gpu<3, halmd::dsfloat> test_suite_type;

            BOOST_DATA_TEST_CASE( ordered, dataset, nparticle ) {
                test_ordered<test_suite_type>(nparticle, nspecies, repeat);
            }
            BOOST_DATA_TEST_CASE( unordered, dataset, nparticle ) {
                test_unordered<test_suite_type>(nparticle, nspecies, repeat);
            }
        BOOST_AUTO_TEST_SUITE_END()
# endif
    BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif
