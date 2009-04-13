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
 * \file   text_ostream_backend.hpp
 * \author Andrey Semashev
 * \date   22.04.2007
 *
 * The header contains implementation of a text output stream sink backend.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_SINKS_TEXT_OSTREAM_BACKEND_HPP_INCLUDED_
#define BOOST_LOG_SINKS_TEXT_OSTREAM_BACKEND_HPP_INCLUDED_

#include <ostream>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/log/detail/prologue.hpp>
#include <boost/log/sinks/basic_sink_backend.hpp>
#include <boost/log/utility/record_writer.hpp>

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

/*!
 * \brief An implementation of a text output stream logging sink backend
 *
 * The sink backend puts formatted log records to one or more text streams.
 * The sink supports \c record_writer interface.
 */
template< typename CharT >
class BOOST_LOG_EXPORT basic_text_ostream_backend :
    public basic_formatting_sink_backend< CharT >
{
    //! Base type
    typedef basic_formatting_sink_backend< CharT > base_type;

public:
    //! Character type
    typedef typename base_type::char_type char_type;
    //! String type to be used as a message text holder
    typedef typename base_type::string_type string_type;
    //! String type to be used as a message text holder
    typedef typename base_type::target_string_type target_string_type;
    //! Attribute values view type
    typedef typename base_type::values_view_type values_view_type;
    //! Log record type
    typedef typename base_type::record_type record_type;
    //! Output stream type
    typedef typename base_type::stream_type stream_type;

private:
    //! \cond

    //! Structure with data regarding a single stream
    struct stream_info
    {
        shared_ptr< stream_type > strm;
        record_writer* record_listener;
    };
    //! Type of the container that holds all aggregated streams
    typedef std::vector< stream_info > ostream_sequence;

    //! \endcond

private:
    //! Output stream list
    ostream_sequence m_Streams;
    //! Auto-flush flag
    bool m_fAutoFlush;

public:
    /*!
     * Constructor. No streams attached to the constructed backend, auto flush feature disabled.
     */
    basic_text_ostream_backend();
    /*!
     * Destructor
     */
    ~basic_text_ostream_backend();

    /*!
     * The method adds a new stream to the sink.
     *
     * \param strm Pointer to the stream. Must not be NULL.
     */
    void add_stream(shared_ptr< stream_type > const& strm);
    /*!
     * The method removes a stream from the sink. If the stream is not attached to the sink,
     * the method has no effect.
     *
     * \param strm Pointer to the stream. Must not be NULL.
     */
    void remove_stream(shared_ptr< stream_type > const& strm);

    /*!
     * Sets the flag to automatically flush buffers of all attached streams after each log record
     */
    void auto_flush(bool f = true);

private:
#ifndef BOOST_LOG_DOXYGEN_PASS
    //! The method writes the message to the sink
    void do_consume(record_type const& record, target_string_type const& formatted_message);
#endif // BOOST_LOG_DOXYGEN_PASS
};

#ifdef BOOST_LOG_USE_CHAR
typedef basic_text_ostream_backend< char > text_ostream_backend;        //!< Convenience typedef for narrow-character logging
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
typedef basic_text_ostream_backend< wchar_t > wtext_ostream_backend;    //!< Convenience typedef for wide-character logging
#endif

} // namespace sinks

} // namespace log

} // namespace boost

#ifdef _MSC_VER
#pragma warning(pop)
#endif // _MSC_VER

#endif // BOOST_LOG_SINKS_TEXT_OSTREAM_BACKEND_HPP_INCLUDED_
