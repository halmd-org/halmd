/*!
 * (C) 2009 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   windows_version.hpp
 * \author Andrey Semashev
 * \date   07.03.2009
 * 
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#ifndef BOOST_LOG_WINDOWS_VERSION_HPP_INCLUDED_
#define BOOST_LOG_WINDOWS_VERSION_HPP_INCLUDED_

#include <boost/log/detail/prologue.hpp>

#ifdef BOOST_WINDOWS

#if defined(BOOST_LOG_USE_WINNT6_API)

#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0600 // _WIN32_WINNT_LONGHORN
#endif

#else

#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0500 // _WIN32_WINNT_WIN2K
#endif

#endif // BOOST_LOG_USE_WINNT6_API

#endif // BOOST_WINDOWS

#endif // BOOST_LOG_WINDOWS_VERSION_HPP_INCLUDED_
