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
 * \file   severity_logger.hpp
 * \author Andrey Semashev
 * \date   08.03.2007
 *
 * The header contains implementation of a logger with severity level support.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_SOURCES_SEVERITY_LOGGER_HPP_INCLUDED_
#define BOOST_LOG_SOURCES_SEVERITY_LOGGER_HPP_INCLUDED_

#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/mpl/aux_/lambda_support.hpp>
#include <boost/log/detail/prologue.hpp>
#include <boost/log/detail/singleton.hpp>
#include <boost/log/detail/new_shared.hpp>
#if !defined(BOOST_LOG_NO_THREADS)
#include <boost/log/detail/thread_specific.hpp>
#include <boost/log/detail/light_rw_mutex.hpp>
#endif
#include <boost/log/sources/basic_logger.hpp>
#include <boost/log/attributes/attribute.hpp>
#include <boost/log/attributes/basic_attribute_value.hpp>
#include <boost/log/keywords/severity.hpp>

#ifdef _MSC_VER
#pragma warning(push)
// 'm_A' : class 'A' needs to have dll-interface to be used by clients of class 'B'
#pragma warning(disable: 4251)
// non dll-interface class 'A' used as base for dll-interface class 'B'
#pragma warning(disable: 4275)
#endif // _MSC_VER

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace sources {

namespace aux {

    //! A helper traits to get severity attribute name constant in the proper type
    template< typename >
    struct severity_attribute_name;

#ifdef BOOST_LOG_USE_CHAR
    template< >
    struct severity_attribute_name< char >
    {
        static const char* get() { return "Severity"; }
    };
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
    template< >
    struct severity_attribute_name< wchar_t >
    {
        static const wchar_t* get() { return L"Severity"; }
    };
#endif

    //! Severity level storage class
    class severity_level_holder :
        public enable_shared_from_this< severity_level_holder >,
        public boost::log::aux::lazy_singleton< severity_level_holder, shared_ptr< severity_level_holder > >
    {
        friend class boost::log::aux::lazy_singleton< severity_level_holder, shared_ptr< severity_level_holder > >;
        typedef boost::log::aux::lazy_singleton< severity_level_holder, shared_ptr< severity_level_holder > > singleton_base;

    private:
#if !defined(BOOST_LOG_NO_THREADS)
        //! The actual severity level value
        boost::log::aux::thread_specific< int > m_Value;
#else
        //! The actual severity level value
        int m_Value;
#endif

    public:
        ~severity_level_holder();

        //! Returns an instance of the holder
        static BOOST_LOG_EXPORT shared_ptr< severity_level_holder > get();

        //! The method sets the actual level
        void set_value(int level)
        {
            m_Value = level;
        }
        //! The method returns the current level
        int get_value() const
        {
#if !defined(BOOST_LOG_NO_THREADS)
            return m_Value.get();
#else
            return m_Value;
#endif
        }

    private:
        severity_level_holder();
        //! Initializes the singleton instance
        static void init_instance();
    };

    //! Severity level attribute implementation
    template< typename LevelT >
    class severity_level :
        public attribute,
        public attribute_value,
        public enable_shared_from_this< severity_level< LevelT > >
    {
    public:
        //! Stored level type
        typedef LevelT held_type;

    private:
        //! Pointer to the level storage
        shared_ptr< severity_level_holder > m_pHolder;

    public:
        //! Default constructor
        severity_level() : m_pHolder(severity_level_holder::get())
        {
        }

        //! The method returns the actual attribute value. It must not return NULL.
        virtual shared_ptr< attribute_value > get_value()
        {
            return this->shared_from_this();
        }
        //! The method sets the actual level
        void set_value(held_type level)
        {
            m_pHolder->set_value(static_cast< int >(level));
        }

        //! The method dispatches the value to the given object
        virtual bool dispatch(type_dispatcher& dispatcher)
        {
            register type_visitor< held_type >* visitor =
                dispatcher.get_visitor< held_type >();
            if (visitor)
            {
                visitor->visit(static_cast< held_type >(m_pHolder->get_value()));
                return true;
            }
            else
                return false;
        }

        //! The method is called when the attribute value is passed to another thread
        virtual shared_ptr< attribute_value > detach_from_thread()
        {
#if !defined(BOOST_LOG_NO_THREADS)
            return boost::log::aux::new_shared<
                attributes::basic_attribute_value< held_type >
            >(static_cast< held_type >(m_pHolder->get_value()));
#else
            // With multithreading disabled we may safely return this here. This method will not be called anyway.
            return this->shared_from_this();
#endif
        }
    };

} // namespace aux

/*!
 * \brief Logger class with severity level support
 *
 * The logger registers a special attribute with an integral value type on construction.
 * This attribute will provide severity level for each log record being made through the logger.
 * The severity level can be omitted on logging record construction, in which case the default
 * level will be used. The default level can also be customized by passing it to the logger constructor.
 */
template< typename BaseT, typename LevelT = int >
class basic_severity_logger :
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
    //! Threading model being used
    typedef typename base_type::threading_model threading_model;
    //! Log record type
    typedef typename base_type::record_type record_type;

    //! Severity level type
    typedef LevelT severity_level;
    //! Severity attribute type
    typedef aux::severity_level< severity_level > severity_attribute;

private:
    //! Default severity
    severity_level m_DefaultSeverity;
    //! Severity attribute
    shared_ptr< severity_attribute > m_pSeverity;

public:
    /*!
     * Default constructor. The constructed logger will have a severity attribute registered.
     * The default level for log records will be 0.
     */
    basic_severity_logger() :
        base_type(),
        m_DefaultSeverity(static_cast< severity_level >(0)),
        m_pSeverity(boost::log::aux::new_shared< severity_attribute >()) // make_shared doesn't work in 1.36: http://svn.boost.org/trac/boost/ticket/2126
    {
        base_type::add_attribute_unlocked(
            aux::severity_attribute_name< char_type >::get(),
            m_pSeverity);
    }
    /*!
     * Copy constructor
     */
    basic_severity_logger(basic_severity_logger const& that) :
        base_type(static_cast< base_type const& >(that)),
        m_DefaultSeverity(that.m_DefaultSeverity),
        m_pSeverity(that.m_pSeverity)
    {
        base_type::attributes()[aux::severity_attribute_name< char_type >::get()] = m_pSeverity;
    }
    /*!
     * Constructor with named arguments. Allows to setup the default level for log records.
     *
     * \param args A set of named arguments. The following arguments are supported:
     *             \li \c severity - default severity value
     */
    template< typename ArgsT >
    explicit basic_severity_logger(ArgsT const& args) :
        base_type(args),
        m_DefaultSeverity(args[keywords::severity | severity_level()]),
        m_pSeverity(boost::log::aux::new_shared< severity_attribute >()) // make_shared doesn't work in 1.36: http://svn.boost.org/trac/boost/ticket/2126
    {
        base_type::add_attribute_unlocked(
            aux::severity_attribute_name< char_type >::get(),
            m_pSeverity);
    }

    /*!
     * The method opens a new logging record with the default severity
     */
    record_type open_record()
    {
        open_record_lock _(this->threading_base());
        return open_record_unlocked();
    }

    /*!
     * The method opens a new logging record. Record level can be specified as one of the named arguments.
     *
     * \param args A set of named arguments. The following arguments are supported:
     *             \li \c severity - log record severity level
     */
    template< typename ArgsT >
    record_type open_record(ArgsT const& args)
    {
        open_record_lock _(this->threading_base());
        return open_record_unlocked(args);
    }

protected:
    /*!
     * Severity attribute accessor
     */
    shared_ptr< severity_attribute > const& severity() const { return m_pSeverity; }
    /*!
     * Default severity value getter
     */
    severity_level default_severity() const { return m_DefaultSeverity; }

    //! Lock requirement for the open_record_unlocked method
    typedef typename strictest_lock<
        typename base_type::open_record_lock,
        no_lock
    >::type open_record_lock;

    /*!
     * Unlocked \c open_record
     */
    record_type open_record_unlocked()
    {
        m_pSeverity->set_value(m_DefaultSeverity);
        return base_type::open_record_unlocked();
    }
    /*!
     * Unlocked \c open_record
     */
    template< typename ArgsT >
    record_type open_record_unlocked(ArgsT const& args)
    {
        m_pSeverity->set_value(args[keywords::severity | m_DefaultSeverity]);
        return base_type::open_record_unlocked();
    }

    //! Lock requirement for the swap_unlocked method
    typedef typename strictest_lock<
        typename base_type::swap_lock,
#ifndef BOOST_LOG_NO_THREADS
        lock_guard< threading_model >
#else
        no_lock
#endif // !defined(BOOST_LOG_NO_THREADS)
    >::type swap_lock;

    //! Unlocked \c swap
    void swap_unlocked(basic_severity_logger& that)
    {
        base_type::swap_unlocked(static_cast< base_type& >(that));
        std::swap(m_DefaultSeverity, that.m_DefaultSeverity);
        m_pSeverity.swap(that.m_pSeverity);
    }

#ifndef BOOST_LOG_DOXYGEN_PASS
public:
    BOOST_MPL_AUX_LAMBDA_SUPPORT(2, basic_severity_logger, (BaseT, LevelT))
#endif // BOOST_LOG_DOXYGEN_PASS
};

#ifndef BOOST_LOG_DOXYGEN_PASS

#ifdef BOOST_LOG_USE_CHAR

//! Narrow-char logger with severity level support
template< typename LevelT = int >
class severity_logger :
    public basic_composite_logger<
        char,
        severity_logger< LevelT >,
        single_thread_model,
        mpl::vector1< basic_severity_logger< mpl::_1, LevelT > >
    >
{
    BOOST_LOG_FORWARD_LOGGER_CONSTRUCTORS_TEMPLATE(severity_logger)
};

#if !defined(BOOST_LOG_NO_THREADS)

//! Narrow-char thread-safe logger with severity level support
template< typename LevelT = int >
class severity_logger_mt :
    public basic_composite_logger<
        char,
        severity_logger_mt< LevelT >,
        multi_thread_model< boost::log::aux::light_rw_mutex >,
        mpl::vector1< basic_severity_logger< mpl::_1, LevelT > >
    >
{
    BOOST_LOG_FORWARD_LOGGER_CONSTRUCTORS_TEMPLATE(severity_logger_mt)
};

#endif // !defined(BOOST_LOG_NO_THREADS)

#endif

#ifdef BOOST_LOG_USE_WCHAR_T

//! Wide-char logger with severity level support
template< typename LevelT = int >
class wseverity_logger :
    public basic_composite_logger<
        wchar_t,
        wseverity_logger< LevelT >,
        single_thread_model,
        mpl::vector1< basic_severity_logger< mpl::_1, LevelT > >
    >
{
    BOOST_LOG_FORWARD_LOGGER_CONSTRUCTORS_TEMPLATE(wseverity_logger)
};

#if !defined(BOOST_LOG_NO_THREADS)

//! Wide-char thread-safe logger with severity level support
template< typename LevelT = int >
class wseverity_logger_mt :
    public basic_composite_logger<
        wchar_t,
        wseverity_logger_mt< LevelT >,
        multi_thread_model< boost::log::aux::light_rw_mutex >,
        mpl::vector1< basic_severity_logger< mpl::_1, LevelT > >
    >
{
    BOOST_LOG_FORWARD_LOGGER_CONSTRUCTORS_TEMPLATE(wseverity_logger_mt)
};

#endif // !defined(BOOST_LOG_NO_THREADS)

#endif

#else // BOOST_LOG_DOXYGEN_PASS

/*!
 * \brief Narrow-char logger. Functionally equivalent to \c basic_severity_logger.
 *
 * See \c basic_severity_logger class template for a more detailed description
 */
template< typename LevelT = int >
class severity_logger :
    public basic_severity_logger<
        basic_logger< char, severity_logger< LevelT >, single_thread_model >,
        LevelT
    >
{
public:
    /*!
     * Default constructor
     */
    severity_logger();
    /*!
     * Copy constructor
     */
    severity_logger(severity_logger const& that);
    /*!
     * Constructor with named arguments
     */
    template< typename... ArgsT >
    explicit severity_logger(ArgsT... const& args);
    /*!
     * Assignment operator
     */
    severity_logger& operator= (severity_logger const& that)
    /*!
     * Swaps two loggers
     */
    void swap(severity_logger& that);
};

/*!
 * \brief Narrow-char thread-safe logger. Functionally equivalent to \c basic_severity_logger.
 *
 * See \c basic_severity_logger class template for a more detailed description
 */
template< typename LevelT = int >
class severity_logger_mt :
    public basic_severity_logger<
        basic_logger< char, severity_logger_mt< LevelT >, multi_thread_model< shared_mutex > >,
        LevelT
    >
{
public:
    /*!
     * Default constructor
     */
    severity_logger_mt();
    /*!
     * Copy constructor
     */
    severity_logger_mt(severity_logger_mt const& that);
    /*!
     * Constructor with named arguments
     */
    template< typename... ArgsT >
    explicit severity_logger_mt(ArgsT... const& args);
    /*!
     * Assignment operator
     */
    severity_logger_mt& operator= (severity_logger_mt const& that)
    /*!
     * Swaps two loggers
     */
    void swap(severity_logger_mt& that);
};

/*!
 * \brief Wide-char logger. Functionally equivalent to \c basic_severity_logger.
 *
 * See \c basic_severity_logger class template for a more detailed description
 */
template< typename LevelT = int >
class wseverity_logger :
    public basic_severity_logger<
        basic_logger< wchar_t, wseverity_logger< LevelT >, single_thread_model >,
        LevelT
    >
{
public:
    /*!
     * Default constructor
     */
    wseverity_logger();
    /*!
     * Copy constructor
     */
    wseverity_logger(wseverity_logger const& that);
    /*!
     * Constructor with named arguments
     */
    template< typename... ArgsT >
    explicit wseverity_logger(ArgsT... const& args);
    /*!
     * Assignment operator
     */
    wseverity_logger& operator= (wseverity_logger const& that)
    /*!
     * Swaps two loggers
     */
    void swap(wseverity_logger& that);
};

/*!
 * \brief Wide-char thread-safe logger. Functionally equivalent to \c basic_severity_logger.
 *
 * See \c basic_severity_logger class template for a more detailed description
 */
template< typename LevelT = int >
class wseverity_logger_mt :
    public basic_severity_logger<
        basic_logger< wchar_t, wseverity_logger_mt< LevelT >, multi_thread_model< shared_mutex > >,
        LevelT
    >
{
public:
    /*!
     * Default constructor
     */
    wseverity_logger_mt();
    /*!
     * Copy constructor
     */
    wseverity_logger_mt(wseverity_logger_mt const& that);
    /*!
     * Constructor with named arguments
     */
    template< typename... ArgsT >
    explicit wseverity_logger_mt(ArgsT... const& args);
    /*!
     * Assignment operator
     */
    wseverity_logger_mt& operator= (wseverity_logger_mt const& that)
    /*!
     * Swaps two loggers
     */
    void swap(wseverity_logger_mt& that);
};

#endif // BOOST_LOG_DOXYGEN_PASS

} // namespace sources

} // namespace log

} // namespace boost

#ifdef _MSC_VER
#pragma warning(pop)
#endif // _MSC_VER

//! The macro allows to put a record with a specific severity level into log
#define BOOST_LOG_STREAM_SEV(logger, lvl)\
    BOOST_LOG_STREAM_WITH_PARAMS((logger), (::boost::log::keywords::severity = (lvl)))

#ifndef BOOST_LOG_NO_SHORTHAND_NAMES

//! An equivalent to BOOST_LOG_STREAM_SEV(logger, lvl)
#define BOOST_LOG_SEV(logger, lvl) BOOST_LOG_STREAM_SEV(logger, lvl)

#endif // BOOST_LOG_NO_SHORTHAND_NAMES

#endif // BOOST_LOG_SOURCES_SEVERITY_LOGGER_HPP_INCLUDED_
