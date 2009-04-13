/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   sink.cpp
 * \author Andrey Semashev
 * \date   03.11.2007
 * 
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#include <list>
#include <boost/none.hpp>
#include <boost/optional.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/log/sinks/sink.hpp>

#if !defined(BOOST_LOG_NO_THREADS)
#include <boost/thread/thread.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace sinks {

namespace aux {

//! Sink implementation data
template< typename CharT >
class asynchronous_sink_impl< CharT >::implementation
{
private:
    //! Synchronization mutex type
    typedef mutex mutex_type;
    //! Lock type
    typedef lock_guard< mutex_type > scoped_lock;

    //! Pending records queue
    typedef std::list< record_type > enqueued_records;

    //! A simple functor to start internal thread and not involve Boost.Bind
    struct thread_starter;
    friend struct thread_starter;
    struct thread_starter
    {
        explicit thread_starter(implementation* pThis) : m_pThis(pThis) {}
        void operator() () const { m_pThis->output_thread(); }

    private:
        implementation* m_pThis;
    };

private:
    //! The flag shows that the output thread should finish
    volatile bool m_Finishing;
    //! Opaque pointer to the sink backend
    void* m_pBackend;
    //! Pointer to the \c consume callback
    consume_callback_t m_ConsumeCallback;

    //! Records queue waiting to be output
    enqueued_records m_EnqueuedRecords;

    //! Synchronization mutex
    mutex_type m_Mutex;
    //! Sleep condition
    condition_variable m_Condition;

    //! Mutex to support locked_backend_ptr
    shared_backend_lock::mutex_type m_SharedBackendMutex;
    //! Shared lock to support locked_backend_ptr
    optional< shared_backend_lock > m_SharedBackendLock;

    //! Record output thread
    optional< thread > m_Thread;

public:
    //! Constructor
    implementation(void* p, consume_callback_t ccb) :
        m_Finishing(false),
        m_pBackend(p),
        m_ConsumeCallback(ccb)
    {
        m_Thread = boost::in_place(thread_starter(this));
    }
    //! Destructor
    ~implementation()
    {
        if (!!m_Thread)
        {
            {
                scoped_lock lock(m_Mutex);
                m_Finishing = true;
            }
            m_Condition.notify_one();
            m_Thread->join();
            m_Thread = none;
        }
    }

    //! The method puts the record into the queue
    void enqueue_message(record_type& record)
    {
        // Make sure that no references to the thread-specific data is left in attribute values
        record.detach_from_thread();

        // Put the record into the queue
        {
            scoped_lock _(m_Mutex);
            m_EnqueuedRecords.push_back(record);
        }

        m_Condition.notify_one();
    }

    //! Accessor to the shared lock for the locked_backend_ptr support
    optional< shared_backend_lock >& get_shared_backend_lock()
    {
        shared_backend_lock::scoped_lock lock(m_SharedBackendMutex);
        if (!m_SharedBackendLock)
            m_SharedBackendLock = boost::in_place(boost::ref(lock));
        return m_SharedBackendLock;
    }

private:
    //! Output thread routine
    void output_thread()
    {
        enqueued_records ExtractedRecords;

        while (!m_Finishing)
        {
            {
                unique_lock< mutex_type > lock(m_Mutex);
                if (m_EnqueuedRecords.empty())
                    m_Condition.wait(lock);

                ExtractedRecords.splice(ExtractedRecords.end(), m_EnqueuedRecords);
            }

            while (!ExtractedRecords.empty())
            {
                try
                {
                    record_type rec;
                    rec.swap(ExtractedRecords.front());
                    ExtractedRecords.pop_front();
                    shared_backend_lock::scoped_lock lock(m_SharedBackendMutex);
                    m_ConsumeCallback(m_pBackend, rec);
                }
                catch (...)
                {
                    // We do nothing here. There's nothing we can do, actually.
                }
            }
        }
    }
};

//! Constructor
template< typename CharT >
asynchronous_sink_impl< CharT >::asynchronous_sink_impl(void* pBackend, consume_callback_t ccb)
    : pImpl(new implementation(pBackend, ccb))
{
}

//! Destructor
template< typename CharT >
asynchronous_sink_impl< CharT >::~asynchronous_sink_impl()
{
    delete pImpl;
}

//! The method puts the record into the queue
template< typename CharT >
void asynchronous_sink_impl< CharT >::enqueue_message(record_type record)
{
    pImpl->enqueue_message(record);
}

//! Accessor to the shared lock for the locked_backend_ptr support
template< typename CharT >
optional< shared_backend_lock >& asynchronous_sink_impl< CharT >::get_shared_backend_lock() const
{
    return pImpl->get_shared_backend_lock();
}

#ifdef BOOST_LOG_USE_CHAR
template class BOOST_LOG_EXPORT asynchronous_sink_impl< char >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
template class BOOST_LOG_EXPORT asynchronous_sink_impl< wchar_t >;
#endif

} // namespace aux

} // namespace sinks

} // namespace log

} // namespace boost

#endif // !defined(BOOST_LOG_NO_THREADS)
