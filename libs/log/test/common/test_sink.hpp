/*!
 * (C) 2009 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * \file   test_sink.hpp
 * \author Andrey Semashev
 * \date   18.03.2009
 *
 * \brief  This header contains a test sink frontend that is used through various tests.
 */

#ifndef BOOST_LOG_TESTS_TEST_SINK_HPP_INCLUDED_
#define BOOST_LOG_TESTS_TEST_SINK_HPP_INCLUDED_

#include <map>
#include <boost/compatibility/cpp_c_headers/cstddef>
#include <boost/log/sinks/sink.hpp>
#include <boost/log/utility/slim_string.hpp>

//! A sink implementation for testing purpose
template< typename CharT >
struct test_sink :
    public boost::log::sinks::sink< CharT >
{
private:
    typedef boost::log::sinks::sink< CharT > base_type;

public:
    typedef typename base_type::char_type char_type;
    typedef typename base_type::string_type string_type;
    typedef typename base_type::values_view_type values_view_type;
    typedef typename base_type::record_type record_type;
    typedef boost::log::basic_slim_string< char_type > slim_string_type;
    typedef std::map< slim_string_type, std::size_t > attr_counters_map;

public:
    attr_counters_map m_Consumed;
    std::size_t m_RecordCounter;

public:
    test_sink() : m_RecordCounter(0) {}

    void consume(record_type const& record)
    {
        ++m_RecordCounter;
        typename values_view_type::const_iterator
            it = record.attribute_values().begin(),
            end = record.attribute_values().end();
        for (; it != end; ++it)
            ++m_Consumed[it->first];
    }

    void clear()
    {
        m_RecordCounter = 0;
        m_Consumed.clear();
    }
};

#endif // BOOST_LOG_TESTS_TEST_SINK_HPP_INCLUDED_
