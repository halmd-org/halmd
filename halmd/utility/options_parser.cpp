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

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <fstream>
#include <iostream>

#include <halmd/utility/options_parser.hpp>

using namespace boost;
using namespace std;

namespace halmd
{

/**
 * setup program options description
 */
options_parser::options_parser(po::options_description const& desc)
  : desc_(desc) {}

/**
 * parse command line options
 */
void options_parser::parse_command_line(int argc, char** argv)
{
    using namespace po::command_line_style;
    po::command_line_parser parser(argc, argv);
    parser.options(desc_);

    // pass an empty positional options description to the command line
    // parser to warn the user of unintentional positional options
    po::positional_options_description pd;
    parser.positional(pd);

    // disallow abbreviated options, which break forward compatibility of
    // user's scripts as new options are added and create ambiguities
    parser.style(default_style & ~allow_guessing);

    po::parsed_options parsed(parser.run());
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
 * parse command line and config file options
 *
 * This is a helper function for main().
 */
void options_parser::parse(int argc, char** argv)
{
    parse_command_line(argc, argv);

    if (vm_.count("config")) {
        vector<string> const& config = vm_["config"].as<vector<string> >();
        for_each(
            config.begin()
          , config.end()
          , lambda::bind(&options_parser::parse_config_file, this, lambda::_1)
        );
    }
}

} // namespace halmd
