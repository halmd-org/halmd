/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <boost/algorithm/string.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/casts.hpp>
#include <boost/lambda/lambda.hpp>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>
#include <iostream>

#include <halmd/options.hpp>
#include <halmd/utility/date_time.hpp>
#include <halmd/utility/luabind.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

namespace halmd
{

// define here to avoid program-wide dependency on <halmd/version.h>
static char const* default_output_file_name = PROGRAM_NAME "_%Y%m%d_%H%M%S";

/**
 * setup program options description
 */
options_parser::options_parser(po::options_description const& desc)
  : desc_(desc)
{
    desc_.add_options()
        ("output,o",
         po::value<string>()->default_value(default_output_file_name)->notifier(
             boost::lambda::bind(
                 &format_local_time
               , boost::lambda::ll_const_cast<string&>(boost::lambda::_1)
               , boost::lambda::_1
             )
         ),
         "output file prefix")
        ("config,C", po::value<vector<string> >(),
         "parameter input file")
        ("trajectory,J", po::value<string>(),
         "trajectory input file")
        ("version",
         "output version and exit")
        ("help",
         "display this help and exit")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    using namespace luabind;
    register_any_converter<string>();
}


/**
 * parse command line options
 */
void options_parser::parse_command_line(int argc, char** argv)
{
    po::command_line_parser parser(argc, argv);
    po::parsed_options parsed(parser.options(desc_).run());
    po::store(parsed, vm_);
    po::notify(vm_);
}

/**
 * parse config file options
 *
 * @param file_name path to configuration file
 */
void options_parser::parse_config_file(std::string const& file_name)
{
    ifstream ifs(file_name.c_str());
    if (ifs.fail()) {
        throw runtime_error("could not open parameter file '" + file_name + "'");
    }
    po::parsed_options parsed(po::parse_config_file(ifs, desc_));
    po::store(parsed, vm_);
    po::notify(vm_);
}

/**
 * print options parser error message to stderr
 *
 * @param error exception deriving from std::exception
 */
void options_parser::print_error(std::exception const& error)
{
    cerr << PROGRAM_NAME ": " << error.what() << endl;
    cerr << "Try `" PROGRAM_NAME " --help' for more information." << endl;
}

/**
 * print options help message to stdout
 */
void options_parser::print_help()
{
    cout << "Usage: " PROGRAM_NAME " [OPTION]..." << endl << endl
         << desc_ << endl;
}

/**
 * print version information to stdout
 */
void options_parser::print_version()
{
    cout << PROJECT_NAME " (" PROGRAM_DESC ") " PROGRAM_VERSION << endl << endl
         << PROGRAM_COPYRIGHT << endl
         << "This is free software. "
            "You may redistribute copies of it under the terms of" << endl
         << "the GNU General Public License "
            "<http://www.gnu.org/licenses/gpl.html>." << endl
         << "There is NO WARRANTY, to the extent permitted by law." << endl;
}


/**
 * parse command line and config file options
 *
 * This is a helper function for main().
 */
void options_parser::parse(int argc, char** argv)
{
    try {
        parse_command_line(argc, argv);
    }
    catch (std::exception const& e) {
        print_error(e);
        throw exit_exception(EXIT_FAILURE);
    }

    if (vm_.count("help")) {
        print_help();
        throw exit_exception(EXIT_SUCCESS);
    }
    if (vm_.count("version")) {
        print_version();
        throw exit_exception(EXIT_SUCCESS);
    }

    try {
        if (vm_.count("config")) {
            vector<string> const& config = vm_["config"].as<vector<string> >();
            for_each(
                config.begin()
              , config.end()
              , bind(&options_parser::parse_config_file, this, _1)
            );
        }
    }
    catch (std::exception const& e) {
        print_error(e);
        throw exit_exception(EXIT_FAILURE);
    }
}

/**
 * register Lua C++ wrapper
 */
static __attribute__((constructor)) void register_lua()
{
    using namespace luabind;
    lua_registry::get()->push_back
    ((
        namespace_("halmd_wrapper")
        [
            class_<options_parser>("options_parser")
                .def("parsed", &options_parser::parsed)
        ]
      , namespace_("boost")
        [
            namespace_("program_options")
            [
                class_<po::options_description>("options_description") //< only register class
              , class_<po::variable_value>("variable_value")
                    .def(constructor<>())
                    .def("empty", &po::variable_value::empty)
                    .def("defaulted", &po::variable_value::defaulted)
                    .def("value", (any const& (po::variable_value::*)() const) &po::variable_value::value)
                    //< only return-by-value is supported by Luabind boost::any converter
            ]
        ]
    ));
}

} // namespace halmd
