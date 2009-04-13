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
 * \brief  An example of logging into a rotating text file.
 *         See the library tutorial for expanded comments on this code.
 *         It may also be worthwhile reading the Wiki requirements page:
 *         http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?Boost.Logging
 */

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif // _CRT_SECURE_NO_WARNINGS

// #define BOOST_LOG_DYN_LINK 1

#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/shared_ptr.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <boost/log/core.hpp>
#include <boost/log/attributes/clock.hpp>
#include <boost/log/attributes/counter.hpp>
#include <boost/log/sources/basic_logger.hpp>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sinks/sink.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/utility/rotating_ofstream.hpp>
#include <boost/log/utility/empty_deleter.hpp>
#include <boost/log/formatters/format.hpp>
#include <boost/log/formatters/attr.hpp>
#include <boost/log/formatters/date_time.hpp>
#include <boost/log/formatters/message.hpp>

namespace logging = boost::log;
namespace attrs = boost::log::attributes;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace fmt = boost::log::formatters;
namespace keywords = boost::log::keywords;

using boost::shared_ptr;

enum { LOG_RECORDS_TO_WRITE = 10000 };

int main(int argc, char* argv[])
{
    try
    {
        // Open a rotating text file
        shared_ptr< std::ostream > strm(
            new logging::rotating_ofstream("test_%N.log", keywords::rotation_size = 16384));
        if (!strm->good())
            throw std::runtime_error("Failed to open a text log file");

        // Create a text file sink
        shared_ptr< sinks::synchronous_sink< sinks::text_ostream_backend > > sink(
            new sinks::synchronous_sink< sinks::text_ostream_backend >);

        sink->locked_backend()->add_stream(strm);

        sink->locked_backend()->set_formatter(
            fmt::format("%1%: [%2%] - %3%")
                % fmt::attr< unsigned int >("Line #")
                % fmt::date_time< boost::posix_time::ptime >("TimeStamp")
                % fmt::message()
            );

        // Add it to the core
        logging::core::get()->add_sink(sink);

        // Add some attributes too
        shared_ptr< logging::attribute > attr(new attrs::local_clock);
        logging::core::get()->add_global_attribute("TimeStamp", attr);
        attr.reset(new attrs::counter< unsigned int >);
        logging::core::get()->add_global_attribute("Line #", attr);

        // Do some logging
        src::logger lg;
        for (unsigned int i = 0; i < LOG_RECORDS_TO_WRITE; ++i)
        {
            BOOST_LOG(lg) << "Some log record";
        }

        return 0;
    }
    catch (std::exception& e)
    {
        std::cout << "FAILURE: " << e.what() << std::endl;
        return 1;
    }
}
