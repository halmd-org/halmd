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
 * \file   to_file.hpp
 * \author Andrey Semashev
 * \date   16.05.2008
 * 
 * The header contains implementation of convenience functions for enabling logging to a file.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_UTILITY_INIT_TO_FILE_HPP_INCLUDED_
#define BOOST_LOG_UTILITY_INIT_TO_FILE_HPP_INCLUDED_

#include <iostream>
#include <fstream>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/parameter/parameters.hpp>
#include <boost/log/detail/prologue.hpp>
#include <boost/log/detail/sink_init_helpers.hpp>
#include <boost/log/core.hpp>
#include <boost/log/sinks/sink.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/utility/empty_deleter.hpp>
#include <boost/log/utility/rotating_ofstream.hpp>
#include <boost/log/keywords/file_name.hpp>

//! \cond
#ifndef BOOST_LOG_NO_THREADS
#define BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL sinks::synchronous_sink
#else
#define BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL sinks::unlocked_sink
#endif
//! \endcond

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace aux {

//! The function constructs the sink and adds it to the core
template< typename CharT, typename ArgsT >
shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::basic_text_ostream_backend< CharT >
    >
> init_log_to_file(shared_ptr< std::basic_ostream< CharT > > const& strm, ArgsT const& args)
{
    typedef sinks::basic_text_ostream_backend< CharT > backend_t;
    shared_ptr< backend_t > pBackend = boost::make_shared< backend_t >();
    pBackend->add_stream(strm);

    aux::setup_formatter(*pBackend, args,
        typename is_void< typename parameter::binding< ArgsT, keywords::tag::format, void >::type >::type());

    pBackend->auto_flush(args[keywords::auto_flush | false]);

    shared_ptr< BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL< backend_t > > pSink =
        boost::make_shared< BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL< backend_t > >(pBackend);

    aux::setup_filter(*pSink, args,
        typename is_void< typename parameter::binding< ArgsT, keywords::tag::filter, void >::type >::type());

    basic_core< CharT >::get()->add_sink(pSink);

    return pSink;
}

//! File stream factory for the case when no rotation is needed
template< typename CharT >
struct ofstream_factory
{
    typedef shared_ptr< std::basic_ostream< CharT > > ostream_ptr;

    template< typename ArgsT >
    static ostream_ptr create(ArgsT const& args)
    {
        return create_impl(
            args[keywords::file_name],
            args[keywords::open_mode | (std::ios_base::out | std::ios_base::trunc)]);
    }

private:
    static ostream_ptr create_impl(const char* file_name, std::ios_base::openmode mode)
    {
        return boost::make_shared< std::basic_ofstream< CharT > >(file_name, mode);
    }
    static ostream_ptr create_impl(std::string const& file_name, std::ios_base::openmode mode)
    {
        return boost::make_shared< std::basic_ofstream< CharT > >(file_name.c_str(), mode);
    }
    static ostream_ptr create_impl(filesystem::path const& file_name, std::ios_base::openmode mode)
    {
        return boost::make_shared< filesystem::basic_ofstream< CharT > >(file_name, mode);
    }

#ifndef BOOST_FILESYSTEM_NARROW_ONLY
    static ostream_ptr create_impl(const wchar_t* file_name, std::ios_base::openmode mode)
    {
        return create_impl(filesystem::wpath(file_name), mode);
    }
    static ostream_ptr create_impl(std::wstring const& file_name, std::ios_base::openmode mode)
    {
        return create_impl(filesystem::wpath(file_name), mode);
    }
    static ostream_ptr create_impl(filesystem::wpath const& file_name, std::ios_base::openmode mode)
    {
        return boost::make_shared< filesystem::basic_ofstream< CharT > >(file_name, mode);
    }
#endif // BOOST_FILESYSTEM_NARROW_ONLY
};

//! File stream factory for the case when rotation is required
template< typename CharT >
struct rotating_ofstream_factory
{
    typedef shared_ptr< std::basic_ostream< CharT > > ostream_ptr;

    template< typename ArgsT >
    static ostream_ptr create(ArgsT const& args)
    {
        return ostream_ptr(new basic_rotating_ofstream< CharT >(args[keywords::file_name], args));
    }
};

//! The function creates a file stream and passes it to another overload to create the sink
template< typename CharT, typename ArgsT >
shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::basic_text_ostream_backend< CharT >
    >
> init_log_to_file(ArgsT const& args)
{
    typedef typename mpl::if_<
        mpl::and_<
            is_void< typename parameter::binding< ArgsT, keywords::tag::rotation_size, void >::type >,
            is_void< typename parameter::binding< ArgsT, keywords::tag::rotation_interval, void >::type >
        >,
        ofstream_factory< CharT >,
        rotating_ofstream_factory< CharT >
    >::type ostream_factory;

    return init_log_to_file(ostream_factory::create(args), args);
}

//! The function wraps the argument into a file_name named argument, if needed
template< typename T >
inline T const& wrap_file_name(T const& arg, mpl::true_)
{
    return arg;
}
template< typename T >
inline typename parameter::aux::tag< keywords::tag::file_name, T const >::type
wrap_file_name(T const& arg, mpl::false_)
{
    return keywords::file_name = arg;
}

} // namespace aux

#ifndef BOOST_LOG_DOXYGEN_PASS

template< typename CharT, typename ArgT1 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1)
{
    return aux::init_log_to_file< CharT >(
        aux::wrap_file_name(arg1, typename parameter::aux::is_named_argument< ArgT1 >::type()));
}

template< typename CharT, typename ArgT1, typename ArgT2 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2)
{
    return aux::init_log_to_file< CharT >(
        (aux::wrap_file_name(arg1, typename parameter::aux::is_named_argument< ArgT1 >::type()), arg2));
}

template< typename CharT, typename ArgT1, typename ArgT2, typename ArgT3 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3)
{
    return aux::init_log_to_file< CharT >(
        (aux::wrap_file_name(arg1, typename parameter::aux::is_named_argument< ArgT1 >::type()), arg2, arg3));
}

template< typename CharT, typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4)
{
    return aux::init_log_to_file< CharT >(
        (aux::wrap_file_name(arg1, typename parameter::aux::is_named_argument< ArgT1 >::type()), arg2, arg3, arg4));
}

template< typename CharT, typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4, typename ArgT5 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4, ArgT5 const& arg5)
{
    return aux::init_log_to_file< CharT >(
        (aux::wrap_file_name(arg1, typename parameter::aux::is_named_argument< ArgT1 >::type()), arg2, arg3, arg4, arg5));
}

template< typename CharT, typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4, typename ArgT5, typename ArgT6 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4, ArgT5 const& arg5, ArgT6 const& arg6)
{
    return aux::init_log_to_file< CharT >(
        (aux::wrap_file_name(arg1, typename parameter::aux::is_named_argument< ArgT1 >::type()), arg2, arg3, arg4, arg5, arg6));
}

template< typename CharT, typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4, typename ArgT5, typename ArgT6, typename ArgT7 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4, ArgT5 const& arg5, ArgT6 const& arg6, ArgT7 const& arg7)
{
    return aux::init_log_to_file< CharT >(
        (aux::wrap_file_name(arg1, typename parameter::aux::is_named_argument< ArgT1 >::type()), arg2, arg3, arg4, arg5, arg6, arg7));
}

#ifdef BOOST_LOG_USE_CHAR

// Overloads for narrow-character logging
template< typename ArgT1 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1)
{
    return init_log_to_file< char >(arg1);
}

template< typename ArgT1, typename ArgT2 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2)
{
    return init_log_to_file< char >(arg1, arg2);
}

template< typename ArgT1, typename ArgT2, typename ArgT3 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3)
{
    return init_log_to_file< char >(arg1, arg2, arg3);
}

template< typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4)
{
    return init_log_to_file< char >(arg1, arg2, arg3, arg4);
}

template< typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4, typename ArgT5 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4, ArgT5 const& arg5)
{
    return init_log_to_file< char >(arg1, arg2, arg3, arg4, arg5);
}

template< typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4, typename ArgT5, typename ArgT6 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4, ArgT5 const& arg5, ArgT6 const& arg6)
{
    return init_log_to_file< char >(arg1, arg2, arg3, arg4, arg5, arg6);
}

template< typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4, typename ArgT5, typename ArgT6, typename ArgT7 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4, ArgT5 const& arg5, ArgT6 const& arg6, ArgT7 const& arg7)
{
    return init_log_to_file< char >(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
}

#endif // BOOST_LOG_USE_CHAR

#ifdef BOOST_LOG_USE_WCHAR_T

// Overloads for wide-character logging
template< typename ArgT1 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> winit_log_to_file(ArgT1 const& arg1)
{
    return init_log_to_file< wchar_t >(arg1);
}

template< typename ArgT1, typename ArgT2 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> winit_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2)
{
    return init_log_to_file< wchar_t >(arg1, arg2);
}

template< typename ArgT1, typename ArgT2, typename ArgT3 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> winit_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3)
{
    return init_log_to_file< wchar_t >(arg1, arg2, arg3);
}

template< typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> winit_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4)
{
    return init_log_to_file< wchar_t >(arg1, arg2, arg3, arg4);
}

template< typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4, typename ArgT5 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> winit_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4, ArgT5 const& arg5)
{
    return init_log_to_file< wchar_t >(arg1, arg2, arg3, arg4, arg5);
}

template< typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4, typename ArgT5, typename ArgT6 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> winit_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4, ArgT5 const& arg5, ArgT6 const& arg6)
{
    return init_log_to_file< wchar_t >(arg1, arg2, arg3, arg4, arg5, arg6);
}

template< typename ArgT1, typename ArgT2, typename ArgT3, typename ArgT4, typename ArgT5, typename ArgT6, typename ArgT7 >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> winit_log_to_file(ArgT1 const& arg1, ArgT2 const& arg2, ArgT3 const& arg3, ArgT4 const& arg4, ArgT5 const& arg5, ArgT6 const& arg6, ArgT7 const& arg7)
{
    return init_log_to_file< wchar_t >(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
}

#endif // BOOST_LOG_USE_WCHAR_T


#else // BOOST_LOG_DOXYGEN_PASS

/*!
 * The function initializes the logging library to write logs to a file stream.
 * 
 * \param args A number of named arguments. The following parameters are supported:
 *             \li \c file_name The file name or its pattern. This parameter is mandatory.
 *             \li \c open_mode The mask that describes the open mode for the file. See <tt>std::ios_base::openmode</tt>.
 *             \li \c rotation_size The size of the file at which rotation should occur. See <tt>basic_rotating_ofstream</tt>.
 *             \li \c rotation_interval The time interval between file rotations. See <tt>basic_rotating_ofstream</tt>.
 *             \li \c filter Specifies a filter to install into the sink. May be a string that represents a filter,
 *                           or a filter lambda expression.
 *             \li \c format Specifies a formatter to install into the sink. May be a string that represents a formatter,
 *                           or a formatter lambda expression (either streaming or Boost.Format-like notation).
 *             \li \c auto_flush A boolean flag that shows whether the sink should automaticallu flush the stream
 *                               after each written record.
 * \return Pointer to the constructed sink.
 * File name to write log to. Must point to a zero-terminated sequence of characters,
 *                  must not be NULL.
 * \return Pointer to the constructed sink.
 */
template< typename CharT, typename... ArgsT >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgsT... const& args);

/*!
 * Equivalent to <tt>init_log_to_file< char >(args...);</tt>
 */
template< typename... ArgsT >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> init_log_to_file(ArgsT... const& args);

/*!
 * Equivalent to <tt>init_log_to_file< wchar_t >(args...);</tt>
 */
template< typename... ArgsT >
inline shared_ptr<
    BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL<
        sinks::text_ostream_backend
    >
> winit_log_to_file(ArgsT... const& args);

#endif // BOOST_LOG_DOXYGEN_PASS

} // namespace log

} // namespace boost

#undef BOOST_LOG_FILE_SINK_FRONTEND_INTERNAL

#endif // BOOST_LOG_UTILITY_INIT_TO_FILE_HPP_INCLUDED_
