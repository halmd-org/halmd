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
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/casts.hpp>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE

#include <halmd/io/logger.hpp>
#include <halmd/script.hpp>
#include <halmd/utility/date_time.hpp>
#include <halmd/utility/hostname.hpp>
#include <halmd/utility/options_parser.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace boost::algorithm;
using namespace halmd;
using namespace std;

/**
 * Run HAL’s MD package
 *
 * This function loads the Lua scripting engine, parses program options
 * from the command line and optionally a config file, sets up logging,
 * and runs the Lua simulation script.
 */
int main(int argc, char **argv)
{
    static logger log;

    log.log_to_console(logger::trace); //< facilitate debugging

    try {
        static script script; //< load Lua script engine

        //
        // assemble program options
        //
        options_parser parser;

        script.options(parser);

        parser.add_options()
            ("output,o",
             po::value<string>()->default_value(PROGRAM_NAME "_%Y%m%d_%H%M%S", "")->notifier(
                 lambda::ll_const_cast<string&>(lambda::_1) = lambda::bind(
                     &format_local_time
                   , lambda::_1
                 )
             ),
             "output file prefix")
            ("config,c", po::value<string>(),
             "parameter input file")
            ("trajectory", po::value<string>(),
             "trajectory input file")
            ("verbose,v", po::accum_value<int>()->default_value(logger::warning),
             "increase verbosity")
            ("version",
             "output version and exit")
            ("help,h",
             "display this help and exit")
            ;

        //
        // parse program options from command line and config file
        //
        po::variables_map vm;
        try {
            parser.parse_command_line(argc, argv, vm);

            if (vm.count("config")) {
                parser.parse_config_file(vm["config"].as<string>(), vm);
            }
        }
        catch (po::error const& e) {
            cerr << PROGRAM_NAME ": " << e.what() << endl;
            cerr << "Try `" PROGRAM_NAME " --help' for more information." << endl;
            return EXIT_FAILURE;
        }

        //
        // print version information to stdout
        //
        if (vm.count("version")) {
            cout << PROJECT_NAME " (" PROGRAM_DESC ") " PROGRAM_VERSION << endl << endl
                 << PROGRAM_COPYRIGHT << endl
                 << "This is free software. "
                    "You may redistribute copies of it under the terms of" << endl
                 << "the GNU General Public License "
                    "<http://www.gnu.org/licenses/gpl.html>." << endl
                 << "There is NO WARRANTY, to the extent permitted by law." << endl;
            return EXIT_SUCCESS;
        }

        //
        // print options help message to stdout
        //
        if (vm.count("help")) {
            cout << "Usage: " PROGRAM_NAME " [OPTION]... [[MODULE] [OPTION]...]..."
                 << endl << endl
                 << parser.options() << endl;
            return EXIT_SUCCESS;
        }

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

        script.parsed(vm); //< pass command line options to Lua

        script.run();
    }
    catch (std::exception const& e) {
        LOG_ERROR(e.what());
        LOG_WARNING(PROJECT_NAME " aborted");
        return EXIT_FAILURE;
    }
    catch (...) {
        LOG_ERROR("unknown exception");
        LOG_WARNING(PROJECT_NAME " aborted");
        return EXIT_FAILURE;
    }

    LOG(PROJECT_NAME " exit");
    return EXIT_SUCCESS;
}
