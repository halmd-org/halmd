/*
 * Copyright © 2008-2013 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/config.hpp>

#include <boost/program_options.hpp>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <cstring>

#include <halmd/io/logger.hpp>
#include <halmd/script.hpp>
#include <halmd/utility/program_options.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/gpu/device.hpp>
#include <halmd/version.h>

using namespace halmd;
using namespace std;

namespace po = boost::program_options;

/**
 * Run HAL’s MD package
 *
 * This function parses the program options, loads the Lua interpreter,
 * stores the command-line arguments in the global table 'arg' with
 * program options at indices < 1 and script options at indices >= 1,
 * and loads the user script from a file or from standard input.
 */
int main(int argc, char **argv)
{
    logging::get().open_console(logging::warning);
    try {
        po::options_description desc;
        desc.add_options()
#ifdef HALMD_WITH_GPU
            ("disable-gpu", "disable GPU acceleration")
            ("gpu-device", po::value<int>()->default_value(-1),
             "GPU device to use")
#endif
            ("help", "display this help and exit")
            ("version", "output version information and exit")
            ;

        vector<string> pos;
        po::variables_map vm;
        try {
            using namespace boost::program_options::command_line_style;
            po::command_line_parser parser(argc, argv);
            parser.extra_style_parser(inject_option_terminator());
            parser.style(default_style & ~allow_guessing);
            po::parsed_options parsed = parser.options(desc).run();
            pos = po::collect_unrecognized(parsed.options, po::include_positional);
            po::store(parsed, vm);
            po::notify(vm);
        }
        catch (po::error const& e) {
            cerr << argv[0] << ": " << e.what() << endl
                 << "Try `" << argv[0] << " --help' for more information." << endl;
            return EXIT_FAILURE;
        }

        if (vm.count("help")) {
            cout << "Usage: " << argv[0] << " [options] [--] script [args]" << endl
                 << "   or: " << argv[0] << " [options] [- [args]]" << endl
                 << endl
                 << desc << endl;
            return EXIT_SUCCESS;
        }

        if (vm.count("version")) {
            cout << PROJECT_NAME " (" PROGRAM_DESC ") " PROGRAM_VERSION << endl;
            if (strlen(PROGRAM_VARIANT) > 0) {
                cout << "variant:" PROGRAM_VARIANT << endl;
            }
            cout << "built with " PROGRAM_COMPILER << " for " << PROGRAM_LIBRARIES << endl << endl
                 << "Copyright © " PROGRAM_COPYRIGHT << endl << endl
                 << "This is free software. "
                    "You may redistribute copies of it under the terms of" << endl
                 << "the GNU Lesser General Public License version 3 or later "
                    "<http://www.gnu.org/licenses/gpl.html>." << endl
                 << "There is NO WARRANTY, to the extent permitted by law." << endl;
            return EXIT_SUCCESS;
        }

        script script;
        luaponte::object arg = luaponte::newtable(script.L);
        luaponte::globals(script.L)["arg"] = arg;

        int offset = -argc + 1;
        if (pos.size() > 1) {
            offset += pos.size() - 1;
        }
        for (int i = 0; i < argc; ++i, ++offset) {
            arg[offset] = string(argv[i]);
        }

#ifdef HALMD_WITH_GPU
        if (vm.count("disable-gpu")) {
#endif
            luaponte::object package = luaponte::globals(script.L)["package"]["loaded"];
            package["halmd.utility.device"] = luaponte::newtable(script.L);
#ifdef HALMD_WITH_GPU
        } else {
            device::set(vm["gpu-device"].as<int>());
        }
#endif

        // read script from file if specified, or from stdin
        if (!pos.empty()) {
            script.dofile(pos.front());
            script.run();
        }
        else {
            script.dofile();
        }
    }
    catch (exception const& e) {
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
