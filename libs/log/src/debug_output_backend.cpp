/*!
 * (C) 2008 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * \file   debug_output_backend.cpp
 * \author Andrey Semashev
 * \date   08.11.2008
 *
 * \brief  A logging sink backend that uses debugger output
 */

#include "windows_version.hpp"
#include <windows.h>
#include <string>
#include <boost/log/sinks/debug_output_backend.hpp>

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace sinks {

namespace {

#if defined(BOOST_LOG_USE_CHAR)
    inline void output_debug_string(const char* str)
    {
        OutputDebugStringA(str);
    }
#endif // defined(BOOST_LOG_USE_CHAR)
#if defined(BOOST_LOG_USE_WCHAR_T)
    inline void output_debug_string(const wchar_t* str)
    {
        OutputDebugStringW(str);
    }
#endif // defined(BOOST_LOG_USE_WCHAR_T)

} // namespace

template< typename CharT >
BOOST_LOG_EXPORT bool basic_debug_output_backend< CharT >::debugger_presence_filter::operator() (values_view_type const& values) const
{
    return (IsDebuggerPresent() != FALSE);
}

template< typename CharT >
basic_debug_output_backend< CharT >::basic_debug_output_backend()
{
}

template< typename CharT >
basic_debug_output_backend< CharT >::~basic_debug_output_backend()
{
}

template< typename CharT >
typename basic_debug_output_backend< CharT >::debugger_presence_filter
basic_debug_output_backend< CharT >::get_debugger_presence_filter() const
{
    return debugger_presence_filter();
}

//! The method puts the formatted message to the event log
template< typename CharT >
void basic_debug_output_backend< CharT >::do_consume(record_type const& record, target_string_type const& formatted_message)
{
    output_debug_string(formatted_message.c_str());
}

#ifdef BOOST_LOG_USE_CHAR
template class BOOST_LOG_EXPORT basic_debug_output_backend< char >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
template class BOOST_LOG_EXPORT basic_debug_output_backend< wchar_t >;
#endif

} // namespace sinks

} // namespace log

} // namespace boost
