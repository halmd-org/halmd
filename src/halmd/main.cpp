/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#include <boost/algorithm/string/join.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/options.hpp>
#include <halmd/script.hpp>
#include <halmd/utility/hostname.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace boost::algorithm;
using namespace halmd;
using namespace std;

/**
 * Run HAL’s MD package
 */
int main(int argc, char **argv)
{
    static logger log;

    log.log_to_console(logger::trace); //< facilitate debugging

    try {
        static script script; //< load Lua script engine

        options_parser options(script.options());
        try {
            options.parse(argc, argv);
        }
        catch (exit_exception const& e) {
            return e.code();
        }
        po::variables_map vm(options.parsed());

        script.init(vm); //< pass command line options to Lua

        log.log_to_console(
            static_cast<logger::severity_level>(vm["verbose"].as<int>())
        );
        log.log_to_file(
            static_cast<logger::severity_level>(
                max(vm["verbose"].as<int>(), static_cast<int>(logger::info))
            )
          , vm["output"].as<string>() + ".log"
        );

        LOG(PROJECT_NAME " (" PROGRAM_DESC ") " PROGRAM_VERSION);
        LOG("variant: " << PROGRAM_VARIANT);
#ifndef NDEBUG
        LOG_WARNING("built with enabled debugging");
#endif
#ifdef __DEVICE_EMULATION__
        LOG_WARNING("built with device emulation");
#endif
        LOG("command line: " << join(vector<string>(argv, argv + argc), " "));
        LOG("host name: " << host_name());

        script.run();
    }
    catch (std::exception const&) {
        LOG_WARNING(PROJECT_NAME " aborted");
        return EXIT_FAILURE;
    }

    LOG(PROJECT_NAME " exit");

    return EXIT_SUCCESS;
}
