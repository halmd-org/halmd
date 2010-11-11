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

#ifndef HALMD_OPTIONS_PARSER_HPP
#define HALMD_OPTIONS_PARSER_HPP

#include <lua.hpp>

#include <halmd/utility/program_options/program_options.hpp>

namespace halmd
{

/**
 * Command-line options parser
 */
class options_parser
{
public:
    static void luaopen(lua_State* L);

    options_parser(po::options_description const& desc);
    void parse_command_line(int argc, char** argv);
    void parse_config_file(std::string const& file_name);
    void print_error(std::exception const& error) const;
    void print_help() const;
    void print_version() const;

    //! returns parsed program options
    po::variables_map const& parsed() const
    {
        return vm_;
    }

    void parse(int argc, char** arg);

private:
    po::options_description desc_; //< options description
    po::variables_map vm_; //< options map
};

/**
 * Exit exception
 */
class exit_exception
{
public:
    //! set exit code
    exit_exception(int code) : code_(code) {}

    //! returns exit code
    int code() const throw()
    {
        return code_;
    }

private:
    int code_; //< exit code
};

} // namespace halmd

#endif /* ! HALMD_OPTIONS_PARSER_HPP */
