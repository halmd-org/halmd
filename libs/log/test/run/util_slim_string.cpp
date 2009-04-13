/*!
 * (C) 2009 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * \file   util_slim_string.cpp
 * \author Andrey Semashev
 * \date   18.01.2009
 *
 * \brief  This header contains tests for the slim_string component.
 */

#define BOOST_TEST_MODULE util_slim_string

#include <cwchar>
#include <set>
#include <vector>
#include <memory>
#include <string>
#include <sstream>
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <boost/log/utility/slim_string.hpp>
#include "char_definitions.hpp"

namespace logging = boost::log;

namespace {

    template< typename CharT >
    inline bool eq_pos(std::size_t ss_pos, std::size_t s_pos)
    {
        if (ss_pos == logging::basic_slim_string< CharT >::npos)
            return (s_pos == std::basic_string< CharT >::npos);
        else
            return (ss_pos == s_pos);
    }

} // namespace

// Tests for constructors
BOOST_AUTO_TEST_CASE_TEMPLATE(constructors, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef typename string_t::traits_type traits_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;

    const CharT* Str = data_t::abc();
    std::size_t Len = traits_t::length(Str);
    const CharT* pStr = Str;

    slim_string_t sstr1 = slim_string_t(Str);
    BOOST_CHECK(traits_t::compare(sstr1.c_str(), Str, Len) == 0);
    BOOST_CHECK_EQUAL(sstr1.length(), Len);
    BOOST_CHECK_EQUAL(sstr1.size(), Len);

    slim_string_t sstr2 = sstr1;
    BOOST_CHECK(traits_t::compare(sstr2.c_str(), Str, Len) == 0);
    BOOST_CHECK_EQUAL(sstr2.length(), Len);
    BOOST_CHECK_EQUAL(sstr2.size(), Len);

    slim_string_t sstr3(pStr, Len - 1);
    BOOST_CHECK(traits_t::compare(sstr3.c_str(), Str, Len - 1) == 0);
    BOOST_CHECK_EQUAL(sstr3.length(), Len - 1);
    BOOST_CHECK_EQUAL(sstr3.size(), Len - 1);

    slim_string_t sstr4 = slim_string_t(string_t(pStr));
    BOOST_CHECK(traits_t::compare(sstr4.c_str(), Str, Len) == 0);
    BOOST_CHECK_EQUAL(sstr4.length(), Len);
    BOOST_CHECK_EQUAL(sstr4.size(), Len);
}

// Tests for copy method
BOOST_AUTO_TEST_CASE_TEMPLATE(copies, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef typename string_t::traits_type traits_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;

    const CharT* pOriginal = data_t::some_test_string();
    std::size_t Len = traits_t::length(pOriginal);
    slim_string_t test = slim_string_t(pOriginal);
    std::vector< CharT > Result(Len + 1, static_cast< CharT >(0));
    test.copy(&Result[0], slim_string_t::npos);
    BOOST_CHECK(traits_t::compare(test.c_str(), &Result[0], Len) == 0);
    BOOST_CHECK_EQUAL(Result[Len], static_cast< CharT >(0));

    std::fill(Result.begin(), Result.end(), static_cast< CharT >(0));
    test.copy(&Result[0], slim_string_t::npos, 5);
    BOOST_CHECK(traits_t::compare(test.c_str() + 5, &Result[0], Len - 5) == 0);
    BOOST_CHECK_EQUAL(Result[Len - 5], static_cast< CharT >(0));

    std::fill(Result.begin(), Result.end(), static_cast< CharT >(0));
    test.copy(&Result[0], 4, 5);
    BOOST_CHECK(traits_t::compare(test.c_str() + 5, &Result[0], 4) == 0);
    BOOST_CHECK_EQUAL(Result[4], static_cast< CharT >(0));
}

// Tests for swapping
BOOST_AUTO_TEST_CASE_TEMPLATE(swapping, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef typename string_t::traits_type traits_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;

    slim_string_t test1 = data_t::abc();
    std::size_t abc_len = traits_t::length(data_t::abc());
    slim_string_t test2 = data_t::some_test_string();
    std::size_t test_string_len = traits_t::length(data_t::some_test_string());

    using std::swap;
    swap(test1, test2);
    BOOST_CHECK(traits_t::compare(test1.c_str(), data_t::some_test_string(), test_string_len) == 0);
    BOOST_CHECK_EQUAL(test1.length(), test_string_len);
    BOOST_CHECK_EQUAL(test1.size(), test_string_len);
    BOOST_CHECK(traits_t::compare(test2.c_str(), data_t::abc(), abc_len) == 0);
    BOOST_CHECK_EQUAL(test2.length(), abc_len);
    BOOST_CHECK_EQUAL(test2.size(), abc_len);

    test1.swap(test2);
    BOOST_CHECK(traits_t::compare(test1.c_str(), data_t::abc(), abc_len) == 0);
    BOOST_CHECK_EQUAL(test1.length(), abc_len);
    BOOST_CHECK_EQUAL(test1.size(), abc_len);
    BOOST_CHECK(traits_t::compare(test2.c_str(), data_t::some_test_string(), test_string_len) == 0);
    BOOST_CHECK_EQUAL(test2.length(), test_string_len);
    BOOST_CHECK_EQUAL(test2.size(), test_string_len);
}

// Tests for streaming output
BOOST_AUTO_TEST_CASE_TEMPLATE(streaming_output, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef typename string_t::traits_type traits_t;
    typedef std::basic_ostringstream< CharT > ostream_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;

    slim_string_t test1 = data_t::some_test_string();

    ostream_t strm;
    strm << test1;
    string_t str = strm.str();
    BOOST_CHECK(str == data_t::some_test_string());
}

// Tests for substr
BOOST_AUTO_TEST_CASE_TEMPLATE(substrs, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef typename string_t::traits_type traits_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;

    slim_string_t test1 = data_t::some_test_string();

    slim_string_t test2 = test1.substr();
    BOOST_CHECK(traits_t::compare(test1.c_str(), test2.c_str(), test1.length()) == 0);
    BOOST_CHECK_EQUAL(test1.length(), test2.length());
    BOOST_CHECK_EQUAL(test1.size(), test2.size());

    slim_string_t test3 = test1.substr(5);
    BOOST_CHECK(traits_t::compare(test1.c_str() + 5, test3.c_str(), test1.length() - 5) == 0);
    BOOST_CHECK_EQUAL(test1.length() - 5, test3.length());
    BOOST_CHECK_EQUAL(test1.size() - 5, test3.size());

    slim_string_t test4 = test1.substr(2, 2);
    BOOST_CHECK(traits_t::compare(test1.c_str() + 2, test4.c_str(), 2) == 0);
    BOOST_CHECK_EQUAL(test4.length(), 2UL);
    BOOST_CHECK_EQUAL(test4.size(), 2UL);
}

// Tests for comparison and ordering
BOOST_AUTO_TEST_CASE_TEMPLATE(comparison, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef typename string_t::traits_type traits_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;

    slim_string_t test = data_t::some_test_string();

    string_t str1 = data_t::some_test_string();
    BOOST_CHECK(test.compare(str1) == 0);

    string_t str2 = str1.substr(5, 4); // results in "test"
    BOOST_CHECK(test.compare(5, 4, str2) == 0);

    BOOST_CHECK(test.compare(str2) < 0);

    string_t str3 = data_t::zero_to_five();
    BOOST_CHECK(test.compare(str3) > 0);

    slim_string_t str4(string_t(data_t::abc() + str1)); // results in "abcsome test string"
    BOOST_CHECK(test.compare(0, slim_string_t::npos, str4, 3, slim_string_t::npos) == 0);

    slim_string_t str5 = data_t::abc();
    BOOST_CHECK(test != str5);

    slim_string_t test1 = test;
    BOOST_CHECK(test == test1);

    slim_string_t test2 = data_t::some_test_string();
    BOOST_CHECK(test == test2);

    // Ordering tests
    std::set< slim_string_t > sstr_set;
    sstr_set.insert(data_t::abc());
    sstr_set.insert(data_t::def());
    sstr_set.insert(data_t::aaa());
    sstr_set.insert(data_t::abcd());
    sstr_set.insert(data_t::zz());
    sstr_set.insert(data_t::ABC());
    sstr_set.insert(data_t::zero_to_five());

    std::set< string_t > str_set;
    str_set.insert(data_t::abc());
    str_set.insert(data_t::def());
    str_set.insert(data_t::aaa());
    str_set.insert(data_t::abcd());
    str_set.insert(data_t::zz());
    str_set.insert(data_t::ABC());
    str_set.insert(data_t::zero_to_five());

    BOOST_CHECK(std::equal(sstr_set.begin(), sstr_set.end(), str_set.begin()));
}

// Tests for find methods
BOOST_AUTO_TEST_CASE_TEMPLATE(finds, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;
    typedef bool (*pos_eq_t)(std::size_t, std::size_t);
    pos_eq_t eq = &eq_pos< CharT >;

    slim_string_t test = data_t::abcdefg0123456789();
    string_t str = data_t::abcdefg0123456789();

    BOOST_CHECK(eq(test.find(data_t::def()), str.find(data_t::def())));
    BOOST_CHECK(eq(test.find(slim_string_t(data_t::zero_to_five())), str.find(data_t::zero_to_five())));
    BOOST_CHECK(eq(test.find(string_t(data_t::abc())), str.find(data_t::abc())));
    BOOST_CHECK(eq(test.find(data_t::some_test_string()), str.find(data_t::some_test_string())));

    BOOST_CHECK(eq(test.find(data_t::abc(), 3), str.find(data_t::abc(), 3)));
    BOOST_CHECK(eq(test.find(data_t::def(), 3), str.find(data_t::def(), 3)));

    BOOST_CHECK(eq(test.find(*str.rbegin()), str.find(*str.rbegin())));
    BOOST_CHECK(eq(test.find(*data_t::zero_to_five(), 5), str.find(*data_t::zero_to_five(), 5)));
}

// Tests for rfind methods
BOOST_AUTO_TEST_CASE_TEMPLATE(rfinds, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;
    typedef bool (*pos_eq_t)(std::size_t, std::size_t);
    pos_eq_t eq = &eq_pos< CharT >;

    slim_string_t test = data_t::abcdefg0123456789();
    string_t str = data_t::abcdefg0123456789();

    BOOST_CHECK(eq(test.rfind(data_t::def()), str.rfind(data_t::def())));
    BOOST_CHECK(eq(test.rfind(slim_string_t(data_t::zero_to_five())), str.rfind(data_t::zero_to_five())));
    BOOST_CHECK(eq(test.rfind(string_t(data_t::abc())), str.rfind(data_t::abc())));
    BOOST_CHECK(eq(test.rfind(data_t::some_test_string()), str.rfind(data_t::some_test_string())));

    BOOST_CHECK(eq(test.rfind(data_t::abc(), 3), str.rfind(data_t::abc(), 3)));
    BOOST_CHECK(eq(test.rfind(data_t::def(), 3), str.rfind(data_t::def(), 3)));

    BOOST_CHECK(eq(test.rfind(*str.rbegin()), str.rfind(*str.rbegin())));
    BOOST_CHECK(eq(test.rfind(*data_t::zero_to_five(), 10), str.rfind(*data_t::zero_to_five(), 10)));
}

// Tests for find_first_of methods
BOOST_AUTO_TEST_CASE_TEMPLATE(find_first_ofs, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;
    typedef bool (*pos_eq_t)(std::size_t, std::size_t);
    pos_eq_t eq = &eq_pos< CharT >;

    slim_string_t test = data_t::abcdefg0123456789();
    string_t str = data_t::abcdefg0123456789();
    string_t a09;
    a09.push_back(str[0]);
    a09.push_back(str[7]);
    a09.push_back(str[16]);
    string_t efg(&str[4], 3);
    string_t fed = data_t::def();
    std::reverse(fed.begin(), fed.end());

    BOOST_CHECK(eq(test.find_first_of(data_t::def()), str.find_first_of(data_t::def())));
    BOOST_CHECK(eq(test.find_first_of(slim_string_t(a09)), str.find_first_of(a09)));
    BOOST_CHECK(eq(test.find_first_of(efg), str.find_first_of(efg)));
    BOOST_CHECK(eq(test.find_first_of(data_t::some_test_string()), str.find_first_of(data_t::some_test_string())));

    BOOST_CHECK(eq(test.find_first_of(data_t::abc(), 5), str.find_first_of(data_t::abc(), 5)));
    BOOST_CHECK(eq(test.find_first_of(fed, 3), str.find_first_of(fed, 3)));

    BOOST_CHECK(eq(test.find_first_of(*str.rbegin()), str.find_first_of(*str.rbegin())));
    BOOST_CHECK(eq(test.find_first_of(*data_t::zero_to_five(), 4), str.find_first_of(*data_t::zero_to_five(), 4)));
}

// Tests for find_last_of methods
BOOST_AUTO_TEST_CASE_TEMPLATE(find_last_ofs, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;
    typedef bool (*pos_eq_t)(std::size_t, std::size_t);
    pos_eq_t eq = &eq_pos< CharT >;

    slim_string_t test = data_t::abcdefg0123456789();
    string_t str = data_t::abcdefg0123456789();
    string_t a09;
    a09.push_back(str[0]);
    a09.push_back(str[7]);
    a09.push_back(str[16]);
    string_t efg(&str[4], 3);
    string_t fed = data_t::def();
    std::reverse(fed.begin(), fed.end());

    BOOST_CHECK(eq(test.find_last_of(data_t::def()), str.find_last_of(data_t::def())));
    BOOST_CHECK(eq(test.find_last_of(slim_string_t(a09)), str.find_last_of(a09)));
    BOOST_CHECK(eq(test.find_last_of(efg), str.find_last_of(efg)));
    BOOST_CHECK(eq(test.find_last_of(data_t::some_test_string()), str.find_last_of(data_t::some_test_string())));

    BOOST_CHECK(eq(test.find_last_of(data_t::abc(), 5), str.find_last_of(data_t::abc(), 5)));
    BOOST_CHECK(eq(test.find_last_of(fed, 3), str.find_last_of(fed, 3)));

    BOOST_CHECK(eq(test.find_last_of(*str.rbegin()), str.find_last_of(*str.rbegin())));
    BOOST_CHECK(eq(test.find_last_of(*data_t::zero_to_five(), 4), str.find_last_of(*data_t::zero_to_five(), 4)));
    BOOST_CHECK(eq(test.find_last_of(str[2], 4), str.find_last_of(str[2], 4)));
}

// Tests for find_first_not_of methods
BOOST_AUTO_TEST_CASE_TEMPLATE(find_first_not_ofs, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;
    typedef bool (*pos_eq_t)(std::size_t, std::size_t);
    pos_eq_t eq = &eq_pos< CharT >;

    slim_string_t test = data_t::abcdefg0123456789();
    string_t str = data_t::abcdefg0123456789();
    string_t ba;
    ba.push_back(str[1]);
    ba.push_back(str[0]);
    string_t a09;
    a09.push_back(str[0]);
    a09.push_back(str[7]);
    a09.push_back(str[16]);
    string_t efg(&str[4], 3);
    string_t fed = data_t::def();
    std::reverse(fed.begin(), fed.end());

    BOOST_CHECK(eq(test.find_first_not_of(ba.c_str()), str.find_first_not_of(ba.c_str())));
    BOOST_CHECK(eq(test.find_first_not_of(slim_string_t(a09)), str.find_first_not_of(a09)));
    BOOST_CHECK(eq(test.find_first_not_of(efg), str.find_first_not_of(efg)));
    BOOST_CHECK(eq(test.find_first_not_of(data_t::some_test_string()), str.find_first_not_of(data_t::some_test_string())));

    BOOST_CHECK(eq(test.find_first_not_of(data_t::abc(), 5), str.find_first_not_of(data_t::abc(), 5)));
    BOOST_CHECK(eq(test.find_first_not_of(fed, 3), str.find_first_not_of(fed, 3)));

    BOOST_CHECK(eq(test.find_first_not_of(*str.rbegin()), str.find_first_not_of(*str.rbegin())));
    BOOST_CHECK(eq(test.find_first_not_of(*data_t::zero_to_five(), 4), str.find_first_not_of(*data_t::zero_to_five(), 4)));
}

// Tests for find_last_not_of methods
BOOST_AUTO_TEST_CASE_TEMPLATE(find_last_not_ofs, CharT, char_types)
{
    typedef test_data< CharT > data_t;
    typedef std::basic_string< CharT > string_t;
    typedef logging::basic_slim_string< CharT > slim_string_t;
    typedef bool (*pos_eq_t)(std::size_t, std::size_t);
    pos_eq_t eq = &eq_pos< CharT >;

    slim_string_t test = data_t::abcdefg0123456789();
    string_t str = data_t::abcdefg0123456789();
    string_t ba98;
    ba98.push_back(str[1]);
    ba98.push_back(str[0]);
    ba98.push_back(str[16]);
    ba98.push_back(str[15]);
    string_t a09;
    a09.push_back(str[0]);
    a09.push_back(str[7]);
    a09.push_back(str[16]);
    string_t efg(&str[4], 3);
    string_t fed = data_t::def();
    std::reverse(fed.begin(), fed.end());

    BOOST_CHECK(eq(test.find_last_not_of(ba98.c_str()), str.find_last_not_of(ba98.c_str())));
    BOOST_CHECK(eq(test.find_last_not_of(slim_string_t(a09)), str.find_last_not_of(a09)));
    BOOST_CHECK(eq(test.find_last_not_of(efg), str.find_last_not_of(efg)));
    BOOST_CHECK(eq(test.find_last_not_of(data_t::some_test_string()), str.find_last_not_of(data_t::some_test_string())));

    BOOST_CHECK(eq(test.find_last_not_of(data_t::abc(), 5), str.find_last_not_of(data_t::abc(), 5)));
    BOOST_CHECK(eq(test.find_last_not_of(fed, 3), str.find_last_not_of(fed, 3)));

    BOOST_CHECK(eq(test.find_last_not_of(*str.rbegin()), str.find_last_not_of(*str.rbegin())));
    BOOST_CHECK(eq(test.find_last_not_of(*data_t::zero_to_five(), 4), str.find_last_not_of(*data_t::zero_to_five(), 4)));
    BOOST_CHECK(eq(test.find_last_not_of(str[2], 4), str.find_last_not_of(str[2], 4)));
}
