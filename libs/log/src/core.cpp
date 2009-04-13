/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * \file   core.cpp
 * \author Andrey Semashev
 * \date   19.04.2007
 *
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#include <memory>
#include <deque>
#include <vector>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/compatibility/cpp_c_headers/cstddef>
#include <boost/log/core.hpp>
#include <boost/log/sinks/sink.hpp>
#include <boost/log/attributes/attribute_values_view.hpp>
#include <boost/log/detail/singleton.hpp>
#if !defined(BOOST_LOG_NO_THREADS)
#include <boost/thread/tss.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/log/detail/shared_lock_guard.hpp>
#include <boost/log/detail/light_rw_mutex.hpp>
#endif

namespace boost {

namespace BOOST_LOG_NAMESPACE {

//! Private record data information, with core-specific structures
template< typename CharT >
struct basic_record< CharT >::private_data :
    public public_data
{
    //! Sink interface type
    typedef sinks::sink< CharT > sink_type;
    //! Sinks container type
    typedef std::deque< shared_ptr< sink_type > > sink_list;

    //! A list of sinks that will accept the record
    sink_list m_AcceptingSinks;

    explicit private_data(values_view_type const& values) : public_data(values)
    {
    }
};

//! Logging system implementation
template< typename CharT >
struct basic_core< CharT >::implementation :
    public log::aux::lazy_singleton<
        implementation,
        shared_ptr< basic_core< CharT > >
    >
{
public:
    //! Base type of singleton holder
    typedef log::aux::lazy_singleton<
        implementation,
        shared_ptr< basic_core< CharT > >
    > base_type;

    //! Front-end class type
    typedef basic_core< char_type > core_type;
#if !defined(BOOST_LOG_NO_THREADS)
    //! Read lock type
    typedef log::aux::shared_lock_guard< log::aux::light_rw_mutex > scoped_read_lock;
    //! Write lock type
    typedef lock_guard< log::aux::light_rw_mutex > scoped_write_lock;
#endif

    //! Sinks container type
    typedef std::vector< shared_ptr< sink_type > > sink_list;

    //! Thread-specific data
    struct thread_data
    {
        //! Thread-specific attribute set
        attribute_set_type ThreadAttributes;
    };

public:
#if !defined(BOOST_LOG_NO_THREADS)
    //! Synchronization mutex
    log::aux::light_rw_mutex Mutex;
#endif

    //! List of sinks involved into output
    sink_list Sinks;

    //! Global attribute set
    attribute_set_type GlobalAttributes;
#if !defined(BOOST_LOG_NO_THREADS)
    //! Thread-specific data
    thread_specific_ptr< thread_data > pThreadData;
#else
    //! Thread-specific data
    std::auto_ptr< thread_data > pThreadData;
#endif

    //! The global state of logging
    volatile bool Enabled;
    //! Global filter
    filter_type Filter;

public:
    //! Constructor
    implementation() : Enabled(true) {}

    //! The method returns the current thread-specific data
    thread_data* get_thread_data()
    {
        thread_data* p = pThreadData.get();
        if (!p)
        {
            init_thread_data();
            p = pThreadData.get();
        }
        return p;
    }

    //! The function initializes the logging system
    static void init_instance()
    {
        base_type::get_instance().reset(new core_type());
    }

private:
    //! The method initializes thread-specific data
    void init_thread_data()
    {
#if !defined(BOOST_LOG_NO_THREADS)
        scoped_write_lock lock(Mutex);
#endif
        if (!pThreadData.get())
        {
            std::auto_ptr< thread_data > p(new thread_data());
            pThreadData.reset(p.get());
            p.release();
        }
    }
};


//! Logging system constructor
template< typename CharT >
basic_core< CharT >::basic_core()
    : pImpl(new implementation())
{
}

//! Logging system destructor
template< typename CharT >
basic_core< CharT >::~basic_core()
{
    delete pImpl;
}

//! The method returns a pointer to the logging system instance
template< typename CharT >
shared_ptr< basic_core< CharT > > basic_core< CharT >::get()
{
    return implementation::get();
}

//! The method enables or disables logging and returns the previous state of logging flag
template< typename CharT >
bool basic_core< CharT >::set_logging_enabled(bool enabled)
{
#if !defined(BOOST_LOG_NO_THREADS)
    typename implementation::scoped_write_lock lock(pImpl->Mutex);
#endif
    const bool old_value = pImpl->Enabled;
    pImpl->Enabled = enabled;
    return old_value;
}

//! The method adds a new sink
template< typename CharT >
void basic_core< CharT >::add_sink(shared_ptr< sink_type > const& s)
{
#if !defined(BOOST_LOG_NO_THREADS)
    typename implementation::scoped_write_lock lock(pImpl->Mutex);
#endif
    typename implementation::sink_list::iterator it =
        std::find(pImpl->Sinks.begin(), pImpl->Sinks.end(), s);
    if (it == pImpl->Sinks.end())
        pImpl->Sinks.push_back(s);
}

//! The method removes the sink from the output
template< typename CharT >
void basic_core< CharT >::remove_sink(shared_ptr< sink_type > const& s)
{
#if !defined(BOOST_LOG_NO_THREADS)
    typename implementation::scoped_write_lock lock(pImpl->Mutex);
#endif
    typename implementation::sink_list::iterator it =
        std::find(pImpl->Sinks.begin(), pImpl->Sinks.end(), s);
    if (it != pImpl->Sinks.end())
        pImpl->Sinks.erase(it);
}


//! The method adds an attribute to the global attribute set
template< typename CharT >
std::pair< typename basic_core< CharT >::attribute_set_type::iterator, bool >
basic_core< CharT >::add_global_attribute(string_type const& name, shared_ptr< attribute > const& attr)
{
#if !defined(BOOST_LOG_NO_THREADS)
    typename implementation::scoped_write_lock lock(pImpl->Mutex);
#endif
    return pImpl->GlobalAttributes.insert(typename attribute_set_type::key_type(name), attr);
}

//! The method removes an attribute from the global attribute set
template< typename CharT >
void basic_core< CharT >::remove_global_attribute(typename attribute_set_type::iterator it)
{
#if !defined(BOOST_LOG_NO_THREADS)
    typename implementation::scoped_write_lock lock(pImpl->Mutex);
#endif
    pImpl->GlobalAttributes.erase(it);
}

//! The method returns the complete set of currently registered global attributes
template< typename CharT >
typename basic_core< CharT >::attribute_set_type basic_core< CharT >::get_global_attributes() const
{
#if !defined(BOOST_LOG_NO_THREADS)
    typename implementation::scoped_write_lock lock(pImpl->Mutex);
#endif
    return pImpl->GlobalAttributes;
}
//! The method replaces the complete set of currently registered global attributes with the provided set
template< typename CharT >
void basic_core< CharT >::set_global_attributes(attribute_set_type const& attrs) const
{
#if !defined(BOOST_LOG_NO_THREADS)
    typename implementation::scoped_write_lock lock(pImpl->Mutex);
#endif
    pImpl->GlobalAttributes = attrs;
}

//! The method adds an attribute to the thread-specific attribute set
template< typename CharT >
std::pair< typename basic_core< CharT >::attribute_set_type::iterator, bool >
basic_core< CharT >::add_thread_attribute(string_type const& name, shared_ptr< attribute > const& attr)
{
    typename implementation::thread_data* p = pImpl->get_thread_data();
    return p->ThreadAttributes.insert(typename attribute_set_type::key_type(name), attr);
}

//! The method removes an attribute from the thread-specific attribute set
template< typename CharT >
void basic_core< CharT >::remove_thread_attribute(typename attribute_set_type::iterator it)
{
    typename implementation::thread_data* p = pImpl->pThreadData.get();
    if (p)
        p->ThreadAttributes.erase(it);
}

//! The method returns the complete set of currently registered thread-specific attributes
template< typename CharT >
typename basic_core< CharT >::attribute_set_type basic_core< CharT >::get_thread_attributes() const
{
    typename implementation::thread_data* p = pImpl->get_thread_data();
    return p->ThreadAttributes;
}
//! The method replaces the complete set of currently registered thread-specific attributes with the provided set
template< typename CharT >
void basic_core< CharT >::set_thread_attributes(attribute_set_type const& attrs) const
{
    typename implementation::thread_data* p = pImpl->get_thread_data();
    p->ThreadAttributes = attrs;
}

//! An internal method to set the global filter
template< typename CharT >
void basic_core< CharT >::set_filter(filter_type const& filter)
{
#if !defined(BOOST_LOG_NO_THREADS)
    typename implementation::scoped_write_lock lock(pImpl->Mutex);
#endif
    pImpl->Filter = filter;
}

//! The method removes the global logging filter
template< typename CharT >
void basic_core< CharT >::reset_filter()
{
#if !defined(BOOST_LOG_NO_THREADS)
    typename implementation::scoped_write_lock lock(pImpl->Mutex);
#endif
    pImpl->Filter.clear();
}

//! The method opens a new record to be written and returns true if the record was opened
template< typename CharT >
typename basic_core< CharT >::record_type basic_core< CharT >::open_record(attribute_set_type const& source_attributes)
{
    record_type rec;

    // Try a quick win first
    if (pImpl->Enabled) try
    {
        typename implementation::thread_data* tsd = pImpl->get_thread_data();
#if !defined(BOOST_LOG_NO_THREADS)
        // Lock the core to be safe against any attribute or sink set modifications
        typename implementation::scoped_read_lock lock(pImpl->Mutex);
#endif

        if (pImpl->Enabled && !pImpl->Sinks.empty())
        {
            // Compose a view of attribute values (unfrozen, yet)
            values_view_type attr_values(source_attributes, tsd->ThreadAttributes, pImpl->GlobalAttributes);

            if (pImpl->Filter.empty() || pImpl->Filter(attr_values))
            {
                // The global filter passed, trying the sinks
                typedef typename record_type::private_data record_private_data;
                record_private_data* pData = NULL;
                typename implementation::sink_list::iterator it = pImpl->Sinks.begin(), end = pImpl->Sinks.end();
                for (; it != end; ++it)
                {
                    try
                    {
                        if (it->get()->will_consume(attr_values))
                        {
                            // If at least one sink accepts the record, it's time to create it
                            if (!pData)
                            {
                                attr_values.freeze();
                                rec.m_pData = pData = new record_private_data(attr_values);
                            }
                            pData->m_AcceptingSinks.push_back(*it);
                        }
                    }
                    catch (...)
                    {
                        // Assume that the sink is incapable to receive messages now
                    }
                }
            }
        }
    }
    catch (...)
    {
        // Something has gone wrong. As the library should impose minimum influence
        // on the user's code, we simply mimic here that the record is not needed.
    }

    return rec;
}

//! The method pushes the record
template< typename CharT >
void basic_core< CharT >::push_record(record_type const& rec)
{
    typedef typename record_type::private_data record_private_data;
    record_private_data* pData = static_cast< record_private_data* >(rec.m_pData.get());

    typename record_private_data::sink_list::iterator
        it = pData->m_AcceptingSinks.begin(),
        end = pData->m_AcceptingSinks.end();
    for (; it != end; ++it) try
    {
        (*it)->consume(rec);
    }
    catch (...)
    {
    }
}

//  Explicitly instantiate core implementation
#ifdef BOOST_LOG_USE_CHAR
template class BOOST_LOG_EXPORT basic_core< char >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
template class BOOST_LOG_EXPORT basic_core< wchar_t >;
#endif

} // namespace log

} // namespace boost
