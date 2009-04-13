/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   severity_level.cpp
 * \author Andrey Semashev
 * \date   10.05.2008
 * 
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#include <boost/log/detail/new_shared.hpp>
#include <boost/log/sources/severity_logger.hpp>

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace sources {

namespace aux {

//! Default constructor
inline severity_level_holder::severity_level_holder()
{
}

//! Destructor
severity_level_holder::~severity_level_holder()
{
}

//! Returns an instance of the attribute
BOOST_LOG_EXPORT shared_ptr< severity_level_holder > severity_level_holder::get()
{
    return singleton_base::get();
}

//! Initializes the singleton instance
void severity_level_holder::init_instance()
{
    singleton_base::get_instance().reset(new severity_level_holder());
}

} // namespace aux

} // namespace sources

} // namespace log

} // namespace boost
