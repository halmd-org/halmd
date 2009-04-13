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
 * \file   ostream.hpp
 * \author Andrey Semashev
 * \date   22.04.2007
 * 
 * The header contains implementation of a hook for streaming formatters.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_FORMATTERS_OSTREAM_HPP_INCLUDED_
#define BOOST_LOG_FORMATTERS_OSTREAM_HPP_INCLUDED_

#include <boost/log/detail/prologue.hpp>
#include <boost/log/formatters/basic_formatters.hpp>
#include <boost/log/formatters/chain.hpp>
#include <boost/log/formatters/wrappers.hpp>

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace formatters {

/*!
 * \brief A placeholder class to represent a stream in lambda expressions of formatters
 * 
 * The \c stream_placeholder class template is a hook that allows to construct formatters
 * from streaming lambda expressions.
 */
template< typename CharT >
struct stream_placeholder
{
    /*!
     * Trap operator to begin building the lambda expression
     * 
     * \param fmt Either a formatter or an object to be wrapped into a formatter
     */
    template< typename T >
    typename wrap_if_not_formatter< CharT, T >::type operator<< (T const& fmt) const
    {
        typedef typename wrap_if_not_formatter< CharT, T >::type result_type;
        return result_type(fmt);
    }

    /*!
     * An overload for C-style strings. These need a special treatment in order to keep
     * a copy of the string in the formatter instead of a raw pointer.
     * 
     * \param s A string constant to be wrapped
     */
    fmt_wrapper< CharT, std::basic_string< CharT > > operator<< (const CharT* s) const
    {
        return fmt_wrapper< CharT, std::basic_string< CharT > >(s);
    }

    static const stream_placeholder instance;
};

template< typename CharT >
const stream_placeholder< CharT > stream_placeholder< CharT >::instance = {};

//  Placeholders to begin lambda expressions
namespace {

#ifdef BOOST_LOG_USE_CHAR
    //! A placeholder used to construct lambda expressions of streaming formatters for narrow-character logging
    stream_placeholder< char > const& ostrm = stream_placeholder< char >::instance;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
    //! A placeholder used to construct lambda expressions of streaming formatters for wide-character logging
    stream_placeholder< wchar_t > const& wostrm = stream_placeholder< wchar_t >::instance;
#endif

} // namespace

/*!
 * The operator chains streaming formatters or a formatter and a wrapped object
 * 
 * \param left Left-hand formatter
 * \param right Right-hand operand. May either be a formatter or some another object. In the latter case the object
 *              will be wrapped into a surrogate formatter that will attempt to put the object into the formatting stream.
 * \return The constructed chained formatters
 */
template< typename CharT, typename LeftFmtT, typename RightT >
inline fmt_chain<
    CharT,
    LeftFmtT,
    typename wrap_if_not_formatter< CharT, RightT >::type
> operator<< (basic_formatter< CharT, LeftFmtT > const& left, RightT const& right)
{
    typedef fmt_chain<
        CharT,
        LeftFmtT,
        typename wrap_if_not_formatter< CharT, RightT >::type
    > result_type;
    return result_type(static_cast< LeftFmtT const& >(left), right);
}

/*!
 * The operator chains a formatter and a C-style string. This case needs a special treatment in order to keep
 * a copy of the string in the formatter instead of a raw pointer.
 * 
 * \param left Left-hand formatter
 * \param str Right-hand string. Must be a zero-terminated sequence of characters, must not be NULL.
 * \return The constructed chained formatters
 */
template< typename FmtT, typename CharT >
inline fmt_chain<
    CharT,
    FmtT,
    fmt_wrapper<
        CharT,
        std::basic_string< CharT >
    >
> operator<< (basic_formatter< CharT, FmtT > const& left, const CharT* str)
{
    return fmt_chain<
        CharT,
        FmtT,
        fmt_wrapper<
            CharT,
            std::basic_string< CharT >
        >
    >(static_cast< FmtT const& >(left), fmt_wrapper< CharT, std::basic_string< CharT > >(str));
}

} // namespace formatters

} // namespace log

} // namespace boost

#endif // BOOST_LOG_FORMATTERS_OSTREAM_HPP_INCLUDED_
