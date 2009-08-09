/* Logging
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_LOG_NO_THREADS
# include <boost/log/attributes/timer.hpp>
# include <boost/log/filters/attr.hpp>
# include <boost/log/formatters/attr.hpp>
# include <boost/log/formatters/date_time.hpp>
# include <boost/log/formatters/format.hpp>
# include <boost/log/formatters/message.hpp>
# include <boost/log/utility/init/common_attributes.hpp>
# include <boost/log/utility/init/to_console.hpp>
# include <boost/log/utility/init/to_file.hpp>
#undef BOOST_LOG_NO_THREADS
#include <ljgpu/util/log.hpp>

#define TIMESTAMP_FORMAT "%d-%m-%Y %H:%M:%S.%f"

namespace logging = boost::log;
namespace fmt = boost::log::formatters;
namespace flt = boost::log::filters;
namespace sinks = boost::log::sinks;
namespace attrs = boost::log::attributes;
namespace src = boost::log::sources;
namespace keywords = boost::log::keywords;

namespace ljgpu { namespace log
{

src::severity_logger<severity_level> slg;

/**
 * initialize logging
 */
void init(std::string const& filename, int verbosity)
{
    // log to file
    logging::init_log_to_file
    (
	filename,
        keywords::format = fmt::format("[%1%] %2%")
            % fmt::date_time("TimeStamp", keywords::format = TIMESTAMP_FORMAT)
            % fmt::message()
    );
    // log to console
    severity_level sl = (verbosity > 1 ? debug : (verbosity > 0 ? info : warning));
    logging::init_log_to_console
    (
	std::clog,
	keywords::filter = flt::attr<severity_level>("Severity") >= sl,
        keywords::format = fmt::format("[%1%] %2%")
            % fmt::date_time("TimeStamp", keywords::format = TIMESTAMP_FORMAT)
            % fmt::message()
    );
    // add timestamp and record counter
    logging::add_common_attributes();
}

}} // namespace ljgpu::log
