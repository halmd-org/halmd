/*!
 * (C) 2009 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * \file   util_static_type_dispatcher.cpp
 * \author Andrey Semashev
 * \date   09.01.2009
 *
 * \brief  This header contains tests for the static type dispatcher.
 */

#define BOOST_TEST_MODULE util_static_type_dispatcher

#include <string>
#include <typeinfo>
#include <boost/mpl/vector.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/log/utility/type_dispatch/static_type_dispatcher.hpp>

namespace logging = boost::log;

namespace {

    // A simple attribute value
    template< typename T >
    struct my_value
    {
        T m_Value;

        explicit my_value(T const& value) : m_Value(value) {}

        // The function passes the contained type into the dispatcher
        bool dispatch(logging::type_dispatcher& dispatcher)
        {
            logging::type_visitor< T >* visitor = dispatcher.get_visitor< T >();
            if (visitor)
            {
                visitor->visit(m_Value);
                return true;
            }
            else
                return false;
        }
    };

    // The function tests general functionality of the type dispatcher
    template< typename DispatcherT >
    void test_general_functionality(DispatcherT& disp)
    {
        // These two attributes are supported by the dispatcher
        my_value< std::string > val1("Hello world!");
        disp.set_expected(val1.m_Value);
        BOOST_CHECK(val1.dispatch(disp));

        my_value< double > val2(1.2);
        disp.set_expected(val2.m_Value);
        BOOST_CHECK(val2.dispatch(disp));

        // This one is not
        my_value< float > val3(static_cast< float >(-4.3));
        disp.set_expected();
        BOOST_CHECK(!val3.dispatch(disp));
    }


    // Type dispatcher for the supported types
    struct my_dispatcher :
        public logging::static_type_dispatcher<
            boost::mpl::vector< int, double, std::string >
        >
    {
        enum type_expected
        {
            none_expected,
            int_expected,
            double_expected,
            string_expected
        };

        my_dispatcher() : m_Expected(none_expected), m_Int(0), m_Double(0.0) {}

        void set_expected()
        {
            m_Expected = none_expected;
        }
        void set_expected(int value)
        {
            m_Expected = int_expected;
            m_Int = value;
        }
        void set_expected(double value)
        {
            m_Expected = double_expected;
            m_Double = value;
        }
        void set_expected(std::string const& value)
        {
            m_Expected = string_expected;
            m_String = value;
        }

        // Implement visitation logic for all supported types
        void visit(int const& value)
        {
            BOOST_CHECK_EQUAL(m_Expected, int_expected);
            BOOST_CHECK_EQUAL(m_Int, value);
        }
        void visit(double const& value)
        {
            BOOST_CHECK_EQUAL(m_Expected, double_expected);
            BOOST_CHECK_CLOSE(m_Double, value, 0.001);
        }
        void visit(std::string const& value)
        {
            BOOST_CHECK_EQUAL(m_Expected, string_expected);
            BOOST_CHECK_EQUAL(m_String, value);
        }

    private:
        type_expected m_Expected;
        int m_Int;
        double m_Double;
        std::string m_String;
    };

} // namespace

// The test checks that general functionality works
BOOST_AUTO_TEST_CASE(type_dispatch)
{
    my_dispatcher disp;
    test_general_functionality(disp);
}

namespace {

    struct my_dispatcher_root :
        public logging::type_dispatcher
    {
        my_dispatcher_root() : m_ExpectedType(0) {}

        void set_expected()
        {
            m_ExpectedType = 0;
        }
        template< typename T >
        void set_expected(T const&)
        {
            m_ExpectedType = &typeid(T);
        }

        std::type_info const* m_ExpectedType;
    };

    struct my_visitor
    {
        template< typename T >
        struct apply :
            public logging::type_visitor< T >
        {
            typedef apply< T > type;

            void visit(T const& value)
            {
                my_dispatcher_root* root = dynamic_cast< my_dispatcher_root* >(this);
                BOOST_REQUIRE(root != 0);
                BOOST_REQUIRE(root->m_ExpectedType != 0);
                BOOST_CHECK(typeid(T) == *root->m_ExpectedType);
            }
        };
    };

    typedef logging::static_type_dispatcher<
        boost::mpl::vector< int, double, std::string >,
        my_visitor,
        my_dispatcher_root
    > my_generated_dispatcher;

} // namespace

// The test checks that type visitors can be generated and type dispatching still works
BOOST_AUTO_TEST_CASE(type_visitors_generation)
{
    my_generated_dispatcher disp;
    test_general_functionality(disp);
}
