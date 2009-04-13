/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   throw_exception.hpp
 * \author Andrey Semashev
 * \date   09.01.2009
 * 
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_DETAIL_THROW_EXCEPTION_HPP_INCLUDED_
#define BOOST_LOG_DETAIL_THROW_EXCEPTION_HPP_INCLUDED_

#include <boost/compatibility/cpp_c_headers/cstdlib>
#include <boost/throw_exception.hpp>
#include <boost/log/detail/prologue.hpp>

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace aux {

//! This function is equivalent to boost::throw_exception,
//! except that it doesn't cause warnings in non-void returning functions
template< typename T >
inline void BOOST_LOG_NORETURN throw_exception(T const& ex)
{
    boost::throw_exception(ex);
    BOOST_LOG_ASSUME(false);
    std::abort();
}

} // namespace aux

} // namespace log

} // namespace boost

#endif // BOOST_LOG_DETAIL_THROW_EXCEPTION_HPP_INCLUDED_
