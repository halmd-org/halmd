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
 * \file   sink.hpp
 * \author Andrey Semashev
 * \date   22.04.2007
 *
 * The header contains implementation of sink frontents and a general interface of sinks
 * for the logging core.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_SINKS_SINK_HPP_INCLUDED_
#define BOOST_LOG_SINKS_SINK_HPP_INCLUDED_

#include <string>
#include <boost/ref.hpp>
#include <boost/none.hpp>
#include <boost/optional.hpp>
#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/function/function1.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/log/detail/prologue.hpp>
#include <boost/log/detail/shared_lock_guard.hpp>
#include <boost/log/record.hpp>
#include <boost/log/sinks/threading_models.hpp>
#include <boost/log/attributes/attribute_values_view.hpp>
#if !defined(BOOST_LOG_NO_THREADS)
#include <boost/detail/atomic_count.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <boost/log/detail/light_rw_mutex.hpp>
#endif

#ifdef _MSC_VER
#pragma warning(push)
// 'm_A' : class 'A' needs to have dll-interface to be used by clients of class 'B'
#pragma warning(disable: 4251)
// non dll-interface class 'A' used as base for dll-interface class 'B'
#pragma warning(disable: 4275)
#endif // _MSC_VER

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace sinks {

//! A base class for a logging sink frontend
template< typename CharT >
class BOOST_LOG_NO_VTABLE sink : noncopyable
{
public:
    //! Character type
    typedef CharT char_type;
    //! String type to be used as a message text holder
    typedef std::basic_string< char_type > string_type;
    //! Log record type
    typedef basic_record< char_type > record_type;
    //! Attribute values view type
    typedef basic_attribute_values_view< char_type > values_view_type;
    //! Filter function type
    typedef function1< bool, values_view_type const& > filter_type;

private:
#if !defined(BOOST_LOG_NO_THREADS)
    //! Mutex type
    typedef boost::log::aux::light_rw_mutex filter_mutex_type;
    //! Read lock type
    typedef boost::log::aux::shared_lock_guard< filter_mutex_type > scoped_read_lock;
    //! Write lock type
    typedef lock_guard< filter_mutex_type > scoped_write_lock;

    //! Synchronization mutex
    filter_mutex_type m_FilterMutex;
#endif
    //! Filter
    filter_type m_Filter;

public:
    virtual ~sink() {}

    /*!
     * The method sets sink-specific filter functional object
     */
    template< typename T >
    void set_filter(T const& filter)
    {
#if !defined(BOOST_LOG_NO_THREADS)
        scoped_write_lock _(m_FilterMutex);
#endif
        m_Filter = filter;
    }
    /*!
     * The method removes the sink-specific filter
     */
    void reset_filter()
    {
#if !defined(BOOST_LOG_NO_THREADS)
        scoped_write_lock _(m_FilterMutex);
#endif
        m_Filter.clear();
    }

    /*!
     * The method returns \c true if no filter is set or the attribute values pass the filter
     *
     * \param attributes A set of attribute values of a logging record
     */
    bool will_consume(values_view_type const& attributes)
    {
#if !defined(BOOST_LOG_NO_THREADS)
        scoped_read_lock _(m_FilterMutex);
#endif
        return (m_Filter.empty() || m_Filter(attributes));
    }

    /*!
     * The method puts logging message to the sink
     *
     * \param record Logging record to consume
     */
    virtual void consume(record_type const& record) = 0;
};

/*!
 * \brief Non-blocking logging sink frontend
 *
 * The sink frontend does not perform thread synchronization and
 * simply passes logging records to the sink backend.
 */
template< typename SinkBackendT >
class unlocked_sink :
    public sink< typename SinkBackendT::char_type >
{
    typedef sink< typename SinkBackendT::char_type > base_type;

public:
    //! Sink implementation type
    typedef SinkBackendT sink_backend_type;
    //! \cond
    BOOST_MPL_ASSERT((is_model_supported< typename sink_backend_type::threading_model, backend_synchronization_tag >));
    //! \endcond

    typedef typename base_type::record_type record_type;
    typedef typename base_type::string_type string_type;

    //! A type of pointer to the backend
    typedef shared_ptr< sink_backend_type > locked_backend_ptr;

private:
    //! Pointer to the sink backend implementation
    shared_ptr< sink_backend_type > m_pBackend;

public:
    /*!
     * Default constructor. Constructs the sink backend instance.
     * Requires the backend to be default-constructible.
     */
    unlocked_sink() : m_pBackend(new sink_backend_type()) {}
    /*!
     * Constructor attaches user-constructed backend instance
     *
     * \param backend Pointer to the backend instance. Must not be NULL.
     */
    explicit unlocked_sink(shared_ptr< sink_backend_type > const& backend) : m_pBackend(backend)
    {
        BOOST_ASSERT(!!m_pBackend);
    }

    /*!
     * Locking accessor to the attached backend
     */
    locked_backend_ptr locked_backend() const { return m_pBackend; }

    /*!
     * The method puts logging message to the sink
     *
     * \param record Logging record to consume
     */
    void consume(record_type const& record)
    {
        m_pBackend->consume(record);
    }
};

#if !defined(BOOST_LOG_NO_THREADS)

namespace aux {

    //! Shared lock object to support locking_ptr
    struct shared_backend_lock
    {
        typedef recursive_mutex mutex_type;
        typedef unique_lock< mutex_type > scoped_lock;

        boost::detail::atomic_count m_RefCounter;
        scoped_lock m_Lock;

        shared_backend_lock(scoped_lock& l) : m_RefCounter(0), m_Lock(boost::move(l))
        {
        }
    };

    //! A pointer type that locks the backend until it's destroyed
    template< typename SinkBackendT >
    class locking_ptr
    {
    public:
        //! Pointed type
        typedef SinkBackendT element_type;

    private:
        //! The pointer to the backend
        shared_ptr< element_type > m_pBackend;
        //! Reference to the shared lock of the backend
        optional< shared_backend_lock >& m_Lock;

    public:
        //! Constructor
        locking_ptr(shared_ptr< SinkBackendT > const& p, optional< shared_backend_lock >& l)
            : m_pBackend(p), m_Lock(l)
        {
            ++m_Lock->m_RefCounter;
        }
        //! Copy constructor
        locking_ptr(locking_ptr const& that) : m_pBackend(that.m_pBackend), m_Lock(that.m_Lock)
        {
            ++m_Lock->m_RefCounter;
        }
        //! Destructor
        ~locking_ptr()
        {
            if (--m_Lock->m_RefCounter == 0)
                m_Lock = none;
        }

        //! Indirection
        element_type* operator-> () const { return m_pBackend.get(); }
        //! Dereferencing
        element_type& operator* () const { return *m_pBackend; }

        //! Accessor to the raw pointer
        element_type* get() const { return m_pBackend.get(); }

    private:
        //! Assignment (closed)
        locking_ptr& operator= (locking_ptr const&);
    };

    //! Free raw pointer getter to assist generic programming
    template< typename SinkBackendT >
    inline SinkBackendT* get_pointer(locking_ptr< SinkBackendT > const& p)
    {
        return p.get();
    }

} // namespace aux

/*!
 * \brief Synchronous logging sink frontend
 *
 * The sink frontend serializes threads before passing logging records to the backend
 */
template< typename SinkBackendT >
class synchronous_sink :
    public sink< typename SinkBackendT::char_type >
{
    typedef sink< typename SinkBackendT::char_type > base_type;

    //! Mutex type
    typedef aux::shared_backend_lock::mutex_type mutex_type;
    //! Lock type
    typedef aux::shared_backend_lock::scoped_lock scoped_lock;

public:
    //! Sink implementation type
    typedef SinkBackendT sink_backend_type;
    //! \cond
    BOOST_MPL_ASSERT((is_model_supported< typename sink_backend_type::threading_model, frontend_synchronization_tag >));
    //! \endcond

    typedef typename base_type::record_type record_type;
    typedef typename base_type::string_type string_type;

#ifndef BOOST_LOG_DOXYGEN_PASS

    //! A pointer type that locks the backend until it's destroyed
    typedef aux::locking_ptr< sink_backend_type > locked_backend_ptr;

#else // BOOST_LOG_DOXYGEN_PASS

    //! A pointer type that locks the backend until it's destroyed
    typedef implementation_defined locked_backend_ptr;

#endif // BOOST_LOG_DOXYGEN_PASS

private:
    //! Synchronization mutex
    mutable mutex_type m_Mutex;
    //! Pointer to the sink backend implementation
    shared_ptr< sink_backend_type > m_pBackend;

    //! A temporary storage to allow locked_backend_ptr to work
    mutable optional< aux::shared_backend_lock > m_SharedBackendLock;

public:
    /*!
     * Default constructor. Constructs the sink backend instance.
     * Requires the backend to be default-constructible.
     */
    synchronous_sink() : m_pBackend(new sink_backend_type()) {}
    /*!
     * Constructor attaches user-constructed backend instance
     *
     * \param backend Pointer to the backend instance. Must not be NULL.
     */
    explicit synchronous_sink(shared_ptr< sink_backend_type > const& backend) : m_pBackend(backend)
    {
        BOOST_ASSERT(!!m_pBackend);
    }

    /*!
     * Locking accessor to the attached backend
     */
    locked_backend_ptr locked_backend() const
    {
        scoped_lock lock(m_Mutex);
        if (!m_SharedBackendLock)
            m_SharedBackendLock = boost::in_place(boost::ref(lock));
        return locked_backend_ptr(m_pBackend, m_SharedBackendLock);
    }

    /*!
     * The method puts logging message to the sink
     *
     * \param record Logging record to consume
     */
    void consume(record_type const& record)
    {
        scoped_lock _(m_Mutex);
        m_pBackend->consume(record);
    }
};

namespace aux {

    //! Asynchronous logging sink implementation
    template< typename CharT >
    class BOOST_LOG_EXPORT asynchronous_sink_impl
    {
    public:
        //! Character type
        typedef CharT char_type;
        //! String type to be used as a message text holder
        typedef std::basic_string< char_type > string_type;
        //! Log record type
        typedef basic_record< char_type > record_type;

        //! Callback type to call \c consume in the backend
        typedef void (*consume_callback_t)(void*, record_type const&);

    private:
        class implementation;
        implementation* pImpl;

    public:
        asynchronous_sink_impl(void* pBackend, consume_callback_t wmc);
        ~asynchronous_sink_impl();

        //! The method puts the record into the queue
        void enqueue_message(record_type record);

        //! Accessor to the shared lock for the locked_backend_ptr support
        optional< shared_backend_lock >& get_shared_backend_lock() const;

    private:
        asynchronous_sink_impl(asynchronous_sink_impl const&);
        asynchronous_sink_impl& operator= (asynchronous_sink_impl const&);
    };

} // namespace aux

/*!
 * \brief Asynchronous logging sink frontend
 *
 * The frontend starts a separate thread on construction. All logging records are passed
 * to the backend in this dedicated thread only.
 */
template< typename SinkBackendT >
class asynchronous_sink :
    public sink< typename SinkBackendT::char_type >
{
    typedef sink< typename SinkBackendT::char_type > base_type;

public:
    //! Sink implementation type
    typedef SinkBackendT sink_backend_type;
    //! \cond
    BOOST_MPL_ASSERT((is_model_supported< typename sink_backend_type::threading_model, single_thread_tag >));
    //! \endcond

    typedef typename base_type::char_type char_type;
    typedef typename base_type::record_type record_type;
    typedef typename base_type::string_type string_type;

#ifndef BOOST_LOG_DOXYGEN_PASS

    //! A pointer type that locks the backend until it's destroyed
    typedef aux::locking_ptr< sink_backend_type > locked_backend_ptr;

#else // BOOST_LOG_DOXYGEN_PASS

    //! A pointer type that locks the backend until it's destroyed
    typedef implementation_defined locked_backend_ptr;

#endif // BOOST_LOG_DOXYGEN_PASS

private:
    //! Pointer to the sink backend implementation
    shared_ptr< sink_backend_type > m_pBackend;
    //! Synchronization mutex
    aux::asynchronous_sink_impl< char_type > m_Impl;

public:
    /*!
     * Default constructor. Constructs the sink backend instance.
     * Requires the backend to be default-constructible.
     */
    asynchronous_sink() :
        m_pBackend(new sink_backend_type()),
        m_Impl(m_pBackend.get(), &asynchronous_sink::consume_trampoline)
    {
    }
    /*!
     * Constructor attaches user-constructed backend instance
     *
     * \param backend Pointer to the backend instance. Must not be NULL.
     */
    explicit asynchronous_sink(shared_ptr< sink_backend_type > const& backend) :
        m_pBackend(backend),
        m_Impl(m_pBackend.get(), &asynchronous_sink::consume_trampoline)
    {
        BOOST_ASSERT(!!backend);
    }

    /*!
     * Locking accessor to the attached backend
     */
    locked_backend_ptr locked_backend() const
    {
        return locked_backend_ptr(m_pBackend, m_Impl.get_shared_backend_lock());
    }

    /*!
     * The method puts logging message to the sink
     *
     * \param record Logging record to consume
     */
    void consume(record_type const& record)
    {
        m_Impl.enqueue_message(record);
    }

private:
#ifndef BOOST_LOG_DOXYGEN_PASS
    //! Trampoline function to invoke the backend
    static void consume_trampoline(void* pBackend, record_type const& record)
    {
        sink_backend_type* p = reinterpret_cast< sink_backend_type* >(pBackend);
        p->consume(record);
    }
#endif // BOOST_LOG_DOXYGEN_PASS
};

#endif // !defined(BOOST_LOG_NO_THREADS)

} // namespace sinks

} // namespace log

} // namespace boost

#ifdef _MSC_VER
#pragma warning(pop)
#endif // _MSC_VER

#endif // BOOST_LOG_SINKS_SINK_HPP_INCLUDED_
