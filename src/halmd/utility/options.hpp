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

#ifndef HALMD_UTILITY_OPTIONS_HPP
#define HALMD_UTILITY_OPTIONS_HPP

#include <boost/program_options.hpp>

#include <halmd/utility/detail/options.hpp>

namespace halmd
{
namespace po
{

/** parsed program options type */
typedef boost::program_options::variables_map options;

using boost::program_options::option;
using boost::program_options::options_description;
using boost::program_options::value;

using boost::program_options::required_option;

struct unparsed_options
{
    /** unrecognised program options from command line parser */
    std::vector<po::option> command_line_options;
    /** unrecognised program options from config file parser */
    std::vector<std::vector<po::option> > config_file_options;
};

extern void parse_options(int argc, char** argv, options& vm);
extern void parse_options(unparsed_options& unparsed, options_description const& opt, options& vm);

class options_parser_error
{
public:
    options_parser_error(int status)
      : status_(status)
    {}

    int status() const throw()
    {
        return status_;
    }

private:
    int status_;
};

} // namespace po

} // namespace halmd

#endif /* ! HALMD_UTILITY_OPTIONS_HPP */
