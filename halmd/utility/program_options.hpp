/*
 * Copyright © 2012  Peter Colberg
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

#ifndef HALMD_UTILITY_PROGRAM_OPTIONS_HPP
#define HALMD_UTILITY_PROGRAM_OPTIONS_HPP

#include <boost/program_options/option.hpp>

namespace halmd {

/**
 * This functor is an extra_style_parser for boost::program_options.
 *
 * If a positional argument is encountered, i.e. a string not beginning
 * with '-', a '--' terminator is prepended to the list of arguments.
 *
 * If the argument equals '-', the argument is substituted with '',
 * and a '--' terminator is prepended to the list of arguments.
 */
class inject_option_terminator
{
private:
    typedef std::vector<boost::program_options::option> result_type;

public:
    result_type operator()(std::vector<std::string>& args) const
    {
        std::vector<std::string>::iterator i(args.begin());
        if ((*i).substr(0, 1) != "-" || *i == "-") {
            if (*i == "-") {
                (*i).clear();
            }
            args.insert(i, "--");
        }
        return result_type();
    }
};

} // namespace halmd

#endif /* ! HALMD_UTILITY_PROGRAM_OPTIONS_HPP */
