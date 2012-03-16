/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#include <boost/program_options.hpp>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE

#include <halmd/io/logger.hpp>
#include <halmd/script.hpp>
#include <halmd/utility/program_options.hpp>
#include <halmd/utility/lua/lua.hpp>
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
    try {
        po::options_description desc;
        desc.add_options()
            ("help,h", "display this help and exit")
            ("version", "output version information and exit")
            ;

        vector<string> pos;
        po::variables_map vm;
        try {
            po::command_line_parser parser(argc, argv);
            parser.extra_style_parser(inject_option_terminator());
            po::parsed_options parsed = parser.options(desc).run();
            pos = po::collect_unrecognized(parsed.options, po::include_positional);
            po::store(parsed, vm);
            po::notify(vm);
        }
        catch (po::error const& e) {
            cerr << PROGRAM_NAME ": " << e.what() << endl
                 << "Try `" PROGRAM_NAME " --help' for more information." << endl;
            return EXIT_FAILURE;
        }

        if (vm.count("help")) {
            cout << "Usage: " PROGRAM_NAME " [options] [--] script [args]" << endl
                 << "   or: " PROGRAM_NAME " [options] [- [args]]" << endl
                 << endl
                 << desc << endl;
            return EXIT_SUCCESS;
        }

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

        script script;
        luabind::object arg = luabind::newtable(script.L);
        luabind::globals(script.L)["arg"] = arg;

        int offset = -argc + 1;
        if (pos.size() > 1) {
            offset += pos.size() - 1;
        }
        for (int i = 0; i < argc; ++i, ++offset) {
            arg[offset] = string(argv[i]);
        }

        string filename = "";
        if (!pos.empty()) {
            filename = *pos.begin();
        }
        script.dofile(filename);
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
