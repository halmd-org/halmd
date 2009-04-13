/*
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * This header is the Boost.Log library implementation, see the library documentation
 * at http://www.boost.org/libs/log/doc/log.html.
 */
/*!
 * \file   from_stream.hpp
 * \author Andrey Semashev
 * \date   22.03.2008
 *
 * The header contains definition of facilities that allows to initialize the library from a
 * settings file.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_UTILITY_INIT_FROM_STREAM_HPP_INCLUDED_
#define BOOST_LOG_UTILITY_INIT_FROM_STREAM_HPP_INCLUDED_

#include <iosfwd>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/function/function1.hpp>
#include <boost/log/detail/prologue.hpp>
#include <boost/log/sinks/sink.hpp>

#ifdef _MSC_VER
#pragma warning(push)
// 'm_A' : class 'A' needs to have dll-interface to be used by clients of class 'B'
#pragma warning(disable: 4251)
// non dll-interface class 'A' used as base for dll-interface class 'B'
#pragma warning(disable: 4275)
#endif // _MSC_VER

namespace boost {

namespace BOOST_LOG_NAMESPACE {

/*!
 * \brief The function registers a factory for a custom sink
 *
 * The function registers a factory for a sink. The factory will be called to create sink
 * instance when the parser discovers the specified sink type in the settings file. The
 * factory must accept a map of parameters [parameter name -> parameter value] that it
 * may use to initialize the sink. The factory must return a non-NULL pointer to the
 * constructed sink instance.
 *
 * \param sink_name The custom sink name. Must point to a zero-terminated sequence of characters,
 *                  must not be NULL.
 * \param factory Custom sink factory function
 */
template< typename CharT >
BOOST_LOG_EXPORT void register_sink_factory(
    const CharT* sink_name,
    function1<
        shared_ptr< sinks::sink< CharT > >,
        std::map< std::basic_string< CharT >, std::basic_string< CharT > > const&
    > const& factory);

/*!
 * \brief The function registers a factory for a custom sink
 *
 * The function registers a factory for a sink. The factory will be called to create sink
 * instance when the parser discovers the specified sink type in the settings file. The
 * factory must accept a map of parameters [parameter name -> parameter value] that it
 * may use to initialize the sink. The factory must return a non-NULL pointer to the
 * constructed sink instance.
 *
 * \param sink_name The custom sink name
 * \param factory Custom sink factory function
 */
template< typename CharT, typename TraitsT, typename AllocatorT >
inline void register_sink_factory(
    std::basic_string< CharT, TraitsT, AllocatorT > const& sink_name,
    function1<
        shared_ptr< sinks::sink< CharT > >,
        std::map< std::basic_string< CharT >, std::basic_string< CharT > > const&
    > const& factory)
{
    register_sink_factory(sink_name.c_str(), factory);
}


/*!
 * The function initializes the logging library from a stream containing logging settings
 *
 * \param strm Stream, that provides library settings
 *
 * \b Throws: An <tt>std::exception</tt>-based exception if the read data cannot be interpreted as the library settings
 */
template< typename CharT >
BOOST_LOG_EXPORT void init_from_stream(std::basic_istream< CharT >& strm);

} // namespace log

} // namespace boost

#ifdef _MSC_VER
#pragma warning(pop)
#endif // _MSC_VER

#endif // BOOST_LOG_UTILITY_INIT_FROM_STREAM_HPP_INCLUDED_
