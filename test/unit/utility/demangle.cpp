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

#define BOOST_TEST_MODULE demangle
#include <boost/test/unit_test.hpp>

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector.hpp>

#include <halmd/utility/demangle.hpp>
#include <test/tools/ctest.hpp>

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

using namespace halmd;
using namespace std;

/**
 * test C++ type demangling
 */
BOOST_AUTO_TEST_CASE( test_demangled_name )
{
#ifdef HALMD_USE_DEMANGLING
    BOOST_TEST_MESSAGE( "compiler implements type name demangling" );
    BOOST_CHECK_EQUAL( demangled_name<test_type<0>::type>(),
                       "main_" );
    BOOST_CHECK_EQUAL( demangled_name(typeid(test_type<0>::type)), // alternative syntax
                       "main_" );
    BOOST_CHECK_EQUAL( demangled_name<test_type<1>::type>(),
                       "core<2>" );
    BOOST_CHECK_EQUAL( demangled_name<test_type<2>::type>(),
                       "core<3>" );
    BOOST_CHECK_EQUAL( demangled_name<test_type<3>::type>(),
                       "mdsim::core<3, float>" );
    BOOST_CHECK_EQUAL( demangled_name<test_type<4>::type>(),
                       "mdsim::gpu::force<3, float>" );
    BOOST_CHECK_EQUAL( demangled_name<test_type<5>::type>(),
                       "mdsim::gpu::force<3, mdsim::core<3, core<3> > >" );
    BOOST_CHECK_EQUAL( demangled_name<test_type<6>::type>(),
                       "mdsim::gpu::force<3, float>::LennardJones<mdsim::gpu::force<3, float> >" );
    BOOST_CHECK_EQUAL( demangled_name<test_type<7>::type>(),
                       "mdsim::gpu::force<3, mdsim::core<3, core<3> > >::LennardJones<mdsim::core<3, core<3> > >" );
    BOOST_CHECK_EQUAL( demangled_name<test_type<8>::type>(),
                       "mdsim::core<3, core<3> >" );
#else /* HALMD_USE_DEMANGLING */
    BOOST_TEST_MESSAGE( "compiler does not implement type name demangling" );
    BOOST_CHECK_EQUAL( demangled_name<test_type<0>::type>(),
                       typeid(test_type<0>::type ).name() );
    BOOST_CHECK_EQUAL( demangled_name(typeid(test_type<0>::type)), // alternative syntax
                       typeid(test_type<0>::type).name() );
    BOOST_CHECK_EQUAL( demangled_name<test_type<1>::type>(),
                       typeid(test_type<1>::type ).name() );
    BOOST_CHECK_EQUAL( demangled_name<test_type<2>::type>(),
                       typeid(test_type<2>::type ).name() );
    BOOST_CHECK_EQUAL( demangled_name<test_type<3>::type>(),
                       typeid(test_type<3>::type ).name() );
    BOOST_CHECK_EQUAL( demangled_name<test_type<4>::type>(),
                       typeid(test_type<4>::type ).name() );
    BOOST_CHECK_EQUAL( demangled_name<test_type<5>::type>(),
                       typeid(test_type<5>::type ).name() );
    BOOST_CHECK_EQUAL( demangled_name<test_type<6>::type>(),
                       typeid(test_type<6>::type ).name() );
    BOOST_CHECK_EQUAL( demangled_name<test_type<7>::type>(),
                       typeid(test_type<7>::type ).name() );
    BOOST_CHECK_EQUAL( demangled_name<test_type<8>::type>(),
                       typeid(test_type<8>::type ).name() );
#endif /* HALMD_USE_DEMANGLING */
}
