/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   main.cpp
 * \author Andrey Semashev
 * \date   16.11.2008
 * 
 * \brief  An example of logging into Windows event log.
 *
 * The example shows the basic usage of the simple Windows NT event log backend.
 * The code defines custom severity levels, initializes the sink and a couple of
 * attributes to test with, and writes several records at different levels.
 * As a result the written records should appear in the Application log, and
 * should be displayed correctly with the Windows event log viewer.
 */

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif // _CRT_SECURE_NO_WARNINGS

#define BOOST_LOG_DYN_LINK 1

#include <stdexcept>
#include <string>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <boost/log/core.hpp>
#include <boost/log/attributes/clock.hpp>
#include <boost/log/attributes/counter.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sinks/sink.hpp>
#include <boost/log/sinks/event_log_backend.hpp>
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

//! Define application-specific severity levels
enum severity_levels
{
    normal,
    warning,
    error
};

int main(int argc, char* argv[])
{
    try
    {
        // Create an event log sink
        shared_ptr< sinks::synchronous_sink< sinks::simple_event_log_backend > > sink(
            new sinks::synchronous_sink< sinks::simple_event_log_backend >());

        sink->locked_backend()->set_formatter(
            fmt::format("%1%: [%2%] - %3%")
                % fmt::attr< unsigned int >("Line #")
                % fmt::date_time< boost::posix_time::ptime >("TimeStamp")
                % fmt::message()
            );

        // We'll have to map our custom levels to the event log event types
        sinks::event_log::custom_event_type_mapping< severity_levels > mapping("Severity");
        mapping[normal] = sinks::event_log::info;
        mapping[warning] = sinks::event_log::warning;
        mapping[error] = sinks::event_log::error;

        sink->locked_backend()->set_event_type_mapper(mapping);

        // Add the sink to the core
        logging::core::get()->add_sink(sink);

        // Add some attributes too
        shared_ptr< logging::attribute > attr(new attrs::local_clock);
        logging::core::get()->add_global_attribute("TimeStamp", attr);
        attr.reset(new attrs::counter< unsigned int >);
        logging::core::get()->add_global_attribute("Line #", attr);

        // Do some logging
        src::severity_logger< severity_levels > lg(keywords::severity = normal);
        BOOST_LOG_SEV(lg, normal) << "Some record for NT event log with normal level";
        BOOST_LOG_SEV(lg, warning) << "Some record for NT event log with warning level";
        BOOST_LOG_SEV(lg, error) << "Some record for NT event log with error level";

        return 0;
    }
    catch (std::exception& e)
    {
        std::cout << "FAILURE: " << e.what() << std::endl;
        return 1;
    }
}
