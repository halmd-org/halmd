/*
 * Copyright © 2010  Peter Colberg
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE demangle
#include <boost/test/unit_test.hpp>

#include <boost/assign/list_of.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/vector.hpp>

#include <halmd/utility/demangle.hpp>

//
// define test templates and types
//

struct main_ {};

template <int dimension>
struct core {};

namespace mdsim {

template <int dimension, typename T>
struct core {};

namespace gpu {

template <int dimension, typename T>
struct force
{
    template <typename U>
    struct LennardJones
    {
        typedef T value_type;
    };
};

}} // namespace mdsim::gpu

template <int i>
struct test_type
{
    typedef typename boost::mpl::at_c<
        boost::mpl::vector<
            main_
          , core<2>
          , core<3>
          , mdsim::core<3, float>
          , mdsim::gpu::force<3, float>
          , mdsim::gpu::force<3, mdsim::core<3, core<3> > >
          , mdsim::gpu::force<3, float>::LennardJones<mdsim::gpu::force<3, float> >
          , mdsim::gpu::force<3, mdsim::core<3, core<3> > >::LennardJones<mdsim::core<3, core<3> > >
          , typename mdsim::gpu::force<3, mdsim::core<3, core<3> > >::LennardJones<mdsim::gpu::force<3, float> >::value_type
        >
      , i
    >::type type;
};

namespace std {

// needed for Boost Test failure message
static inline ostream& operator<<(ostream& os, vector<string> const& v)
{
    for_each(v.begin(), v.end(), os << boost::lambda::_1);
    return os;
}

} // namespace std

using namespace boost::assign;
using namespace halmd;
using namespace std;

/**
 * test C++ type demangling
 */
BOOST_AUTO_TEST_CASE( test_demangled_name )
{
    BOOST_CHECK_EQUAL( demangled_name<test_type<0>::type >(),
                       "main_");
    BOOST_CHECK_EQUAL( demangled_name(typeid(test_type<0>::type)), // alternative syntax
                       "main_");
    BOOST_CHECK_EQUAL( demangled_name<test_type<1>::type >(),
                       "core<2>");
    BOOST_CHECK_EQUAL( demangled_name<test_type<2>::type >(),
                       "core<3>");
    BOOST_CHECK_EQUAL( demangled_name<test_type<3>::type >(),
                       "mdsim::core<3, float>");
    BOOST_CHECK_EQUAL( demangled_name<test_type<4>::type >(),
                       "mdsim::gpu::force<3, float>");
    BOOST_CHECK_EQUAL( demangled_name<test_type<5>::type >(),
                       "mdsim::gpu::force<3, mdsim::core<3, core<3> > >");
    BOOST_CHECK_EQUAL( demangled_name<test_type<6>::type >(),
                       "mdsim::gpu::force<3, float>::LennardJones<mdsim::gpu::force<3, float> >");
    BOOST_CHECK_EQUAL( demangled_name<test_type<7>::type >(),
                       "mdsim::gpu::force<3, mdsim::core<3, core<3> > >::LennardJones<mdsim::core<3, core<3> > >");
    BOOST_CHECK_EQUAL( demangled_name<test_type<8>::type >(),
                       "mdsim::core<3, core<3> >");
}

/**
 * test type tokenization without template parameters
 */
BOOST_AUTO_TEST_CASE( test_tokenized_name )
{
    BOOST_CHECK_EQUAL( tokenized_name<test_type<0>::type >(),
                       list_of("main_") );
    BOOST_CHECK_EQUAL( tokenized_name<test_type<1>::type >(),
                       list_of("core") );
    BOOST_CHECK_EQUAL( tokenized_name<test_type<2>::type >(),
                       list_of("core") );
    BOOST_CHECK_EQUAL( tokenized_name<test_type<3>::type >(),
                       list_of("mdsim")("core") );
    BOOST_CHECK_NE( tokenized_name<test_type<3>::type >(), // validate vector comparison
                       list_of("mdsimcore") );
    BOOST_CHECK_EQUAL( tokenized_name<test_type<4>::type >(),
                       list_of("mdsim")("gpu")("force") );
    BOOST_CHECK_EQUAL( tokenized_name<test_type<5>::type >(),
                       list_of("mdsim")("gpu")("force") );
    BOOST_CHECK_EQUAL( tokenized_name<test_type<6>::type >(),
                       list_of("mdsim")("gpu")("force")("LennardJones") );
    BOOST_CHECK_EQUAL( tokenized_name<test_type<7>::type >(),
                       list_of("mdsim")("gpu")("force")("LennardJones") );
    BOOST_CHECK_EQUAL( tokenized_name(typeid(test_type<7>::type)), // alternative syntax
                       list_of("mdsim")("gpu")("force")("LennardJones") );
    BOOST_CHECK_EQUAL( tokenized_name<test_type<8>::type >(),
                       list_of("mdsim")("core") );
}
