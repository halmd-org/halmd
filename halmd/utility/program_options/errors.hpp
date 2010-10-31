/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_UTILITY_PROGRAM_OPTIONS_ERRORS_HPP
#define HALMD_UTILITY_PROGRAM_OPTIONS_ERRORS_HPP

#include <boost/program_options.hpp>

namespace halmd
{
namespace po
{

class conflicting_option
  : public boost::program_options::error
{
public:
    conflicting_option(std::string const& first, std::string const& second)
      : boost::program_options::error("option " + first + " conflicts with " + second)
      , option_name_(first) {}

    ~conflicting_option() throw() {}

    std::string const& get_option_name() const throw()
    {
        return option_name_;
    }

private:
    std::string option_name_;
};

class dependent_option
  : public boost::program_options::error
{
public:
    dependent_option(std::string const& first, std::string const& second)
      : boost::program_options::error("option " + first + " depends on " + second)
      , option_name_(first) {}

    ~dependent_option() throw() {}

    std::string const& get_option_name() const throw()
    {
        return option_name_;
    }

private:
    std::string option_name_;
};

} // namespace po

} // namespace halmd

#endif /* ! HALMD_UTILITY_PROGRAM_OPTIONS_ERRORS_HPP */
