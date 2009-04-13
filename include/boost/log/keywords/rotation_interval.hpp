/*
 * (C) 2009 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * This header is the Boost.Log library implementation, see the library documentation
 * at http://www.boost.org/libs/log/doc/log.html.
 */
/*!
 * \file   keywords/rotation_interval.hpp
 * \author Andrey Semashev
 * \date   14.03.2009
 *
 * The header contains the \c rotation_interval keyword declaration.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_KEYWORDS_ROTATION_INTERVAL_HPP_INCLUDED_
#define BOOST_LOG_KEYWORDS_ROTATION_INTERVAL_HPP_INCLUDED_

#include <boost/parameter/keyword.hpp>
#include <boost/log/detail/prologue.hpp>

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace keywords {

    //! The keyword allows to pass maximum time interval of using the same log file to the rotating file stream methods
    BOOST_PARAMETER_KEYWORD(tag, rotation_interval)

} // namespace keywords

} // namespace log

} // namespace boost

#endif // BOOST_LOG_KEYWORDS_ROTATION_INTERVAL_HPP_INCLUDED_
