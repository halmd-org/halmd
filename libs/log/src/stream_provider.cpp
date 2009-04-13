/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   stream_provider.cpp
 * \author Andrey Semashev
 * \date   17.04.2008
 * 
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#include <memory>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/detail/singleton.hpp>
#if !defined(BOOST_LOG_NO_THREADS)
#include <boost/thread/tss.hpp>
#endif

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace aux {

namespace {

//! The pool of stream compounds
template< typename CharT >
class stream_compound_pool :
    public log::aux::lazy_singleton<
        stream_compound_pool< CharT >,
#if !defined(BOOST_LOG_NO_THREADS)
        thread_specific_ptr< stream_compound_pool< CharT > >
#else
        std::auto_ptr< stream_compound_pool< CharT > >
#endif
    >
{
    //! Self type
    typedef stream_compound_pool< CharT > this_type;
#if !defined(BOOST_LOG_NO_THREADS)
    //! Thread-specific pointer type
    typedef thread_specific_ptr< this_type > tls_ptr_type;
#else
    //! Thread-specific pointer type
    typedef std::auto_ptr< this_type > tls_ptr_type;
#endif
    //! Singleton base type
    typedef log::aux::lazy_singleton<
        this_type,
        tls_ptr_type
    > base_type;
    //! Stream compound type
    typedef typename stream_provider< CharT >::stream_compound stream_compound_t;

public:
    //! Pooled stream compounds
    stream_compound_t* m_Top;

    ~stream_compound_pool()
    {
        register stream_compound_t* p = NULL;
        while ((p = m_Top) != NULL)
        {
            m_Top = p->next;
            delete p;
        }
    }

    //! The method returns pool instance
    static stream_compound_pool& get()
    {
        tls_ptr_type& ptr = base_type::get();
        register this_type* p = ptr.get();
        if (!p)
        {
            std::auto_ptr< this_type > pNew(new this_type());
            ptr.reset(pNew.get());
            p = pNew.release();
        }
        return *p;
    }

private:
    stream_compound_pool() : m_Top(NULL) {}
};

} // namespace

//! The method returns an allocated stream compound
template< typename CharT >
typename stream_provider< CharT >::stream_compound* stream_provider< CharT >::allocate_compound(record_type const& rec)
{
    stream_compound_pool< char_type >& pool = stream_compound_pool< char_type >::get();
    if (pool.m_Top)
    {
        register stream_compound* p = pool.m_Top;
        pool.m_Top = p->next;
        p->next = NULL;
        p->stream.record(rec);
        return p;
    }
    else
        return new stream_compound(rec);
}

//! The method releases a compound
template< typename CharT >
void stream_provider< CharT >::release_compound(stream_compound* compound) /* throw() */
{
    stream_compound_pool< char_type >& pool = stream_compound_pool< char_type >::get();
    compound->next = pool.m_Top;
    pool.m_Top = compound;
    compound->stream.record(record_type());
}

//! Explicitly instantiate stream_provider implementation
#ifdef BOOST_LOG_USE_CHAR
template struct stream_provider< char >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
template struct stream_provider< wchar_t >;
#endif

} // namespace aux

} // namespace log

} // namespace boost
