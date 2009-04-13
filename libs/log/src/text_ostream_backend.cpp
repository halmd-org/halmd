/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   text_ostream_backend.cpp
 * \author Andrey Semashev
 * \date   19.04.2007
 * 
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#include <algorithm>
#include <boost/log/sinks/text_ostream_backend.hpp>

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace sinks {

namespace {

    //! A simple lookup predicate
    template< typename StreamT >
    struct stream_lookup_fun
    {
        typedef bool return_type;

        explicit stream_lookup_fun(shared_ptr< StreamT > const& strm) : m_Stream(strm) {}

        template< typename T >
        return_type operator() (T const& that) const
        {
            return (that.strm == m_Stream);
        }

    private:
        shared_ptr< StreamT > const& m_Stream;
    };

    template< typename StreamT >
    static inline stream_lookup_fun< StreamT > stream_lookup(shared_ptr< StreamT > const& strm)
    {
        return stream_lookup_fun< StreamT >(strm);
    }

} // namespace

//! Constructor
template< typename CharT >
basic_text_ostream_backend< CharT >::basic_text_ostream_backend() : m_fAutoFlush(false)
{
}

//! Destructor (just to make it link from the shared library)
template< typename CharT >
basic_text_ostream_backend< CharT >::~basic_text_ostream_backend() {}

//! The method adds a new stream to the sink
template< typename CharT >
void basic_text_ostream_backend< CharT >::add_stream(shared_ptr< stream_type > const& strm)
{
    typename ostream_sequence::iterator it = std::find_if(m_Streams.begin(), m_Streams.end(), stream_lookup(strm));
    if (it == m_Streams.end())
    {
        stream_info info = { strm, dynamic_cast< record_writer* >(strm.get()) };
        m_Streams.push_back(info);
    }
}

//! The method removes a stream from the sink
template< typename CharT >
void basic_text_ostream_backend< CharT >::remove_stream(shared_ptr< stream_type > const& strm)
{
    typename ostream_sequence::iterator it = std::find_if(m_Streams.begin(), m_Streams.end(), stream_lookup(strm));
    if (it != m_Streams.end())
        m_Streams.erase(it);
}

//! Sets the flag to automatically flush buffers after each logged line
template< typename CharT >
void basic_text_ostream_backend< CharT >::auto_flush(bool f)
{
    m_fAutoFlush = f;
}

//! The method writes the message to the sink
template< typename CharT >
void basic_text_ostream_backend< CharT >::do_consume(
    record_type const& record, target_string_type const& message)
{
    typename string_type::const_pointer const p = message.data();
    typename string_type::size_type const s = message.size();
    typename ostream_sequence::const_iterator it = m_Streams.begin();
    for (; it != m_Streams.end(); ++it)
    {
        register stream_type* const strm = it->strm.get();
        if (strm->good()) try
        {
            if (it->record_listener)
                it->record_listener->on_start_record();

            strm->write(p, static_cast< std::streamsize >(s));
            (*strm) << '\n';

            if (it->record_listener)
                it->record_listener->on_end_record();

            if (m_fAutoFlush)
                strm->flush();
        }
        catch (std::exception&)
        {
        }
    }
}

//! Explicitly instantiate sink backend implementation
#ifdef BOOST_LOG_USE_CHAR
template class BOOST_LOG_EXPORT basic_text_ostream_backend< char >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
template class BOOST_LOG_EXPORT basic_text_ostream_backend< wchar_t >;
#endif

} // namespace sinks

} // namespace log

} // namespace boost
