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
#define BOOST_TEST_MODULE rank
#include <boost/test/unit_test.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/mpl/range_c.hpp>

#include <halmd/utility/module/rank.hpp>

using namespace boost;
using namespace boost::lambda;
using namespace halmd::utility;
using namespace std;

template <int i>
struct test_module
{};

template <>
struct test_module<4>
  : test_module<7>
{
    typedef test_module<7> _Base;
};

template <>
struct test_module<6>
  : test_module<9>
{
    typedef test_module<9> _Base;
};

template <>
struct test_module<5>
  : test_module<6>
{
    typedef test_module<6> _Base;
};

/**
 * Transform mpl::int_ into module rank
 */
template <typename T>
struct make_rank_type
{
    typedef new_ptr<module::typed_rank<test_module<T::value> > > type;
};

/**
 * Construct vector of module ranks
 */
struct make_rank_vector
{
    make_rank_vector()
    {
        mpl::for_each<mpl::range_c<int, 0, 10>, make_rank_type<mpl::_1> >(ref(*this));
    }

    template <typename T>
    void operator()(T const& new_ptr)
    {
        rank_vector.push_back(module::rank(new_ptr()));
    }

    std::vector<module::rank> rank_vector;
};

BOOST_FIXTURE_TEST_CASE( test_rank_vector, make_rank_vector )
{
    for (size_t i = 0; i < rank_vector.size(); ++i) {
        BOOST_CHECK_EQUAL(rank_vector[i]->name(), "test_module<" + lexical_cast<string>(i) + ">");
    }
}

/**
 * sort elements by (inverse) collation order
 */
template <typename T>
struct compare_collation_order
{
    bool operator()(T const& left, T const& right) const
    {
        if (typeid(*left) != typeid(*right)) {
            return !typeid(*left).before(typeid(*right));
        }
        return false;
    }
};

BOOST_FIXTURE_TEST_CASE( test_rank_set, make_rank_vector )
{
    set<module::rank, module::compare_rank> rank_set;
    set<module::rank, compare_collation_order<module::rank> > rank_collation_set;
    BOOST_FOREACH( module::rank const& rank, rank_vector ) {
        BOOST_TEST_MESSAGE("insert " + rank->name());
        BOOST_CHECK(rank_set.insert(rank).second);
        BOOST_CHECK(!rank_set.insert(rank).second);
        BOOST_CHECK(rank_collation_set.insert(rank).second);
        BOOST_CHECK(!rank_collation_set.insert(rank).second);

        // FIXME compare sets

        BOOST_FOREACH( module::rank const& rank, rank_set ) {
            BOOST_TEST_MESSAGE(rank->name());
        }
        BOOST_FOREACH( module::rank const& rank, rank_collation_set ) {
            BOOST_TEST_MESSAGE(rank->name());
        }
    }
}
