/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   main.cpp
 * \author Andrey Semashev
 * \date   26.04.2008
 * 
 * \brief  An example of initializing the library from a settings file.
 *         See the library tutorial for expanded comments on this code.
 *         It may also be worthwhile reading the Wiki requirements page:
 *         http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?Boost.Logging
 */

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif // _CRT_SECURE_NO_WARNINGS

// #define BOOST_LOG_DYN_LINK 1

#include <exception>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/shared_ptr.hpp>

#include <boost/log/core.hpp>
#include <boost/log/attributes/clock.hpp>
#include <boost/log/attributes/constant.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/init/from_stream.hpp>
#include <boost/log/utility/empty_deleter.hpp>
#include <boost/log/utility/scoped_attribute.hpp>

// Auto-link static libs
#define BOOST_LIB_NAME boost_system
#include <boost/config/auto_link.hpp>
#define BOOST_LIB_NAME boost_filesystem
#include <boost/config/auto_link.hpp>
#define BOOST_LIB_NAME boost_regex
#include <boost/config/auto_link.hpp>
#define BOOST_LIB_NAME boost_date_time
#include <boost/config/auto_link.hpp>

namespace logging = boost::log;
namespace attrs = boost::log::attributes;
namespace src = boost::log::sources;

using boost::shared_ptr;

// Here we define our application severity levels.
enum severity_level
{
    normal,
    notification,
    warning,
    error,
    critical
};

//  Global logger declaration
BOOST_LOG_DECLARE_GLOBAL_LOGGER(test_lg, src::severity_logger_mt< >)

void try_logging()
{
    src::severity_logger_mt< >& lg = get_test_lg();
    BOOST_LOG_SEV(lg, normal) << "This is a normal severity record";
    BOOST_LOG_SEV(lg, notification) << "This is a notification severity record";
    BOOST_LOG_SEV(lg, warning) << "This is a warning severity record";
    BOOST_LOG_SEV(lg, error) << "This is a error severity record";
    BOOST_LOG_SEV(lg, critical) << "This is a critical severity record";
}

int main(int argc, char* argv[])
{
    try
    {
        // Open the file
        std::ifstream settings("settings.txt");
        if (!settings.is_open())
        {
            std::cout << "Could not open settings.txt file" << std::endl;
            return 1;
        }

        // Read the settings and initialize logging library
        logging::init_from_stream(settings);

        // Add some attributes
        shared_ptr< logging::attribute > attr(new attrs::local_clock);
        logging::core::get()->add_global_attribute("TimeStamp", attr);

        // Try logging
        try_logging();

        // Now enable tagging and try again
        BOOST_LOG_SCOPED_THREAD_TAG("Tag", std::string, "TAGGED");
        try_logging();

        return 0;
    }
    catch (std::exception& e)
    {
        std::cout << "FAILURE: " << e.what() << std::endl;
        return 1;
    }
}
