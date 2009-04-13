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
 * \file   record_writer.hpp
 * \author Andrey Semashev
 * \date   04.01.2008
 * 
 * The header contains definition of an additional interface for a stream to receive
 * additional meta-information about log records during the output. This interface
 * is also supported in the text_ostream_backend. 
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_UTILITY_RECORD_WRITER_HPP_INCLUDED_
#define BOOST_LOG_UTILITY_RECORD_WRITER_HPP_INCLUDED_

#include <boost/log/detail/prologue.hpp>

namespace boost {

namespace BOOST_LOG_NAMESPACE {

/*!
 * \brief A simple interface to extend stream's awareness of the data flow
 * 
 * The stream class should derive from this class in order to be able to receive
 * notifications about logging record boundaries in the written data flow.
 * 
 * This interface is used by the \c basic_text_ostream_backend sink backend.
 */
struct BOOST_LOG_EXPORT record_writer
{
    /*!
     * Destructor
     */
    virtual ~record_writer() {}

    /*!
     * The method is called before any data that belong to the record is written to the stream
     */
    virtual void on_start_record() {}
    /*!
     * The method is called after all data of the record is written to the stream
     */
    virtual void on_end_record() {}
};

} // namespace log

} // namespace boost

#endif // BOOST_LOG_UTILITY_RECORD_WRITER_HPP_INCLUDED_
