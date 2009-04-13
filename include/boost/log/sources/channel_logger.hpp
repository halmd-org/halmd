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
 * \file   channel_logger.hpp
 * \author Andrey Semashev
 * \date   28.02.2008
 *
 * The header contains implementation of a logger with channel support.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_SOURCES_CHANNEL_LOGGER_HPP_INCLUDED_
#define BOOST_LOG_SOURCES_CHANNEL_LOGGER_HPP_INCLUDED_

#include <boost/shared_ptr.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/aux_/lambda_support.hpp>
#include <boost/parameter/binding.hpp>
#include <boost/type_traits/is_void.hpp>
#include <boost/log/detail/prologue.hpp>
#include <boost/log/detail/new_shared.hpp>
#include <boost/log/detail/light_rw_mutex.hpp>
#include <boost/log/sources/basic_logger.hpp>
#include <boost/log/keywords/channel.hpp>
#include <boost/log/attributes/constant.hpp>

#ifdef _MSC_VER
#pragma warning(push)
// 'm_A' : class 'A' needs to have dll-interface to be used by clients of class 'B'
#pragma warning(disable: 4251)
#endif // _MSC_VER

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace sources {

namespace aux {

    //! A helper traits to get channel attribute name constant in the proper type
    template< typename >
    struct channel_attribute_name;

#ifdef BOOST_LOG_USE_CHAR
    template< >
    struct channel_attribute_name< char >
    {
        static const char* get() { return "Channel"; }
    };
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
    template< >
    struct channel_attribute_name< wchar_t >
    {
        static const wchar_t* get() { return L"Channel"; }
    };
#endif

} // namespace aux

/*!
 * \brief Logger class with channel support
 *
 * The logger automatically registers constant attribute with the specified on construction
 * string value, which is a channel name. The channel name cannot be modified through the logger
 * life time.
 */
template< typename BaseT, typename ChannelT = typename BaseT::string_type >
class basic_channel_logger :
    public BaseT
{
    //! Base type
    typedef BaseT base_type;

public:
    //! Character type
    typedef typename base_type::char_type char_type;
    //! Final type
    typedef typename base_type::final_type final_type;
    //! Attribute set type
    typedef typename base_type::attribute_set_type attribute_set_type;
    //! String type
    typedef typename base_type::string_type string_type;
    //! Threading model being used
    typedef typename base_type::threading_model threading_model;

    //! Channel type
    typedef ChannelT channel_type;
    //! Channel attribute type
    typedef attributes::constant< channel_type > channel_attribute;

private:
    //! Channel attribute
    shared_ptr< channel_attribute > m_pChannel;

public:
    /*!
     * Default constructor. The constructed logger does not have the channel attribute.
     */
    basic_channel_logger() : base_type()
    {
    }
    /*!
     * Copy constructor
     */
    basic_channel_logger(basic_channel_logger const& that) :
        base_type(static_cast< base_type const& >(that)),
        m_pChannel(that.m_pChannel)
    {
    }
    /*!
     * Constructor with arguments. Allows to register a channel name attribute on construction.
     *
     * \param args A set of named arguments. The following arguments are supported:
     *             \li \c channel - a string that represents the channel name
     */
    template< typename ArgsT >
    explicit basic_channel_logger(ArgsT const& args) :
        base_type(args)
    {
        init_channel_attribute(args, typename is_void<
            typename parameter::binding< ArgsT, keywords::tag::channel, void >::type
        >::type());
    }

protected:
    //! Lock requirement for the swap_unlocked method
    typedef typename strictest_lock<
        typename base_type::swap_lock,
#ifndef BOOST_LOG_NO_THREADS
        lock_guard< threading_model >
#else
        no_lock
#endif // !defined(BOOST_LOG_NO_THREADS)
    >::type swap_lock;

    /*!
     * Unlocked swap
     */
    void swap_unlocked(basic_channel_logger& that)
    {
        base_type::swap_unlocked(static_cast< base_type& >(that));
        m_pChannel.swap(that.m_pChannel);
    }

private:
#ifndef BOOST_LOG_DOXYGEN_PASS
    //! Initializes the channel attribute
    template< typename ArgsT >
    void init_channel_attribute(ArgsT const& args, mpl::false_ const&)
    {
        channel_type channel_name(args[keywords::channel]);
        m_pChannel = boost::log::aux::new_shared< channel_attribute >(channel_name);
        base_type::add_attribute_unlocked(
            aux::channel_attribute_name< char_type >::get(),
            m_pChannel);
    }
    //! Initializes the channel attribute (dummy, if no channel is specified)
    template< typename ArgsT >
    void init_channel_attribute(ArgsT const& args, mpl::true_ const&)
    {
    }

public:
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1, basic_channel_logger, (BaseT))
#endif // BOOST_LOG_DOXYGEN_PASS
};

#ifndef BOOST_LOG_DOXYGEN_PASS

#ifdef BOOST_LOG_USE_CHAR

//! Narrow-char logger with channel support
template< typename ChannelT = std::string >
class channel_logger :
    public basic_composite_logger<
        char,
        channel_logger< ChannelT >,
        single_thread_model,
        mpl::vector1< basic_channel_logger< mpl::_1, ChannelT > >
    >
{
    BOOST_LOG_FORWARD_LOGGER_CONSTRUCTORS_TEMPLATE(channel_logger)
};

#if !defined(BOOST_LOG_NO_THREADS)

//! Narrow-char thread-safe logger with channel support
template< typename ChannelT = std::string >
class channel_logger_mt :
    public basic_composite_logger<
        char,
        channel_logger_mt< ChannelT >,
        multi_thread_model< boost::log::aux::light_rw_mutex >,
        mpl::vector1< basic_channel_logger< mpl::_1, ChannelT > >
    >
{
    BOOST_LOG_FORWARD_LOGGER_CONSTRUCTORS_TEMPLATE(channel_logger_mt)
};

#endif // !defined(BOOST_LOG_NO_THREADS)

#endif // BOOST_LOG_USE_CHAR

#ifdef BOOST_LOG_USE_WCHAR_T

//! Wide-char logger with channel support
template< typename ChannelT = std::wstring >
class wchannel_logger :
    public basic_composite_logger<
        wchar_t,
        wchannel_logger< ChannelT >,
        single_thread_model,
        mpl::vector1< basic_channel_logger< mpl::_1, ChannelT > >
    >
{
    BOOST_LOG_FORWARD_LOGGER_CONSTRUCTORS_TEMPLATE(wchannel_logger)
};

#if !defined(BOOST_LOG_NO_THREADS)

//! Wide-char thread-safe logger with channel support
template< typename ChannelT = std::wstring >
class wchannel_logger_mt :
    public basic_composite_logger<
        wchar_t,
        wchannel_logger< ChannelT >,
        multi_thread_model< boost::log::aux::light_rw_mutex >,
        mpl::vector1< basic_channel_logger< mpl::_1, ChannelT > >
    >
{
    BOOST_LOG_FORWARD_LOGGER_CONSTRUCTORS_TEMPLATE(wchannel_logger_mt)
};

#endif // !defined(BOOST_LOG_NO_THREADS)

#endif // BOOST_LOG_USE_WCHAR_T

#else // BOOST_LOG_DOXYGEN_PASS

/*!
 * \brief Narrow-char logger. Functionally equivalent to \c basic_channel_logger.
 *
 * See \c basic_channel_logger class template for a more detailed description
 */
template< typename ChannelT = std::string >
class channel_logger :
    public basic_channel_logger<
        basic_logger< char, channel_logger< ChannelT >, single_thread_model >,
        ChannelT
    >
{
public:
    /*!
     * Default constructor
     */
    channel_logger();
    /*!
     * Copy constructor
     */
    channel_logger(channel_logger const& that);
    /*!
     * Constructor with named arguments
     */
    template< typename... ArgsT >
    explicit channel_logger(ArgsT... const& args);
    /*!
     * Assignment operator
     */
    channel_logger& operator= (channel_logger const& that)
    /*!
     * Swaps two loggers
     */
    void swap(channel_logger& that);
};

/*!
 * \brief Narrow-char thread-safe logger. Functionally equivalent to \c basic_channel_logger.
 *
 * See \c basic_channel_logger class template for a more detailed description
 */
template< typename ChannelT = std::string >
class channel_logger_mt :
    public basic_channel_logger<
        basic_logger< char, channel_logger_mt< ChannelT >, multi_thread_model< shared_mutex > >,
        ChannelT
    >
{
public:
    /*!
     * Default constructor
     */
    channel_logger_mt();
    /*!
     * Copy constructor
     */
    channel_logger_mt(channel_logger_mt const& that);
    /*!
     * Constructor with named arguments
     */
    template< typename... ArgsT >
    explicit channel_logger_mt(ArgsT... const& args);
    /*!
     * Assignment operator
     */
    channel_logger_mt& operator= (channel_logger_mt const& that)
    /*!
     * Swaps two loggers
     */
    void swap(channel_logger_mt& that);
};

/*!
 * \brief Wide-char logger. Functionally equivalent to \c basic_channel_logger.
 *
 * See \c basic_channel_logger class template for a more detailed description
 */
template< typename ChannelT = std::wstring >
class wchannel_logger :
    public basic_channel_logger<
        basic_logger< wchar_t, wchannel_logger< ChannelT >, single_thread_model >,
        ChannelT
    >
{
public:
    /*!
     * Default constructor
     */
    wchannel_logger();
    /*!
     * Copy constructor
     */
    wchannel_logger(wchannel_logger const& that);
    /*!
     * Constructor with named arguments
     */
    template< typename... ArgsT >
    explicit wchannel_logger(ArgsT... const& args);
    /*!
     * Assignment operator
     */
    wchannel_logger& operator= (wchannel_logger const& that)
    /*!
     * Swaps two loggers
     */
    void swap(wchannel_logger& that);
};

/*!
 * \brief Wide-char thread-safe logger. Functionally equivalent to \c basic_channel_logger.
 *
 * See \c basic_channel_logger class template for a more detailed description
 */
template< typename ChannelT = std::wstring >
class wchannel_logger_mt :
    public basic_channel_logger<
        basic_logger< wchar_t, wchannel_logger_mt< ChannelT >, multi_thread_model< shared_mutex > >,
        ChannelT
    >
{
public:
    /*!
     * Default constructor
     */
    wchannel_logger_mt();
    /*!
     * Copy constructor
     */
    wchannel_logger_mt(wchannel_logger_mt const& that);
    /*!
     * Constructor with named arguments
     */
    template< typename... ArgsT >
    explicit wchannel_logger_mt(ArgsT... const& args);
    /*!
     * Assignment operator
     */
    wchannel_logger_mt& operator= (wchannel_logger_mt const& that)
    /*!
     * Swaps two loggers
     */
    void swap(wchannel_logger_mt& that);
};

#endif // BOOST_LOG_DOXYGEN_PASS

} // namespace sources

} // namespace log

} // namespace boost

#ifdef _MSC_VER
#pragma warning(pop)
#endif // _MSC_VER

#endif // BOOST_LOG_SOURCES_CHANNEL_LOGGER_HPP_INCLUDED_
