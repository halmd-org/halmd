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
    options_parser() {}
    po::options_description_easy_init add_options();
    void add(po::options_description const& desc);
    void add(po::options_description const& desc, std::string const& section);
    po::options_description options() const;

    void parse_command_line(std::vector<std::string> const& args, po::variables_map& vm);
    void parse_command_line(int argc, char** argv, po::variables_map& vm);
    void parse_config_file(std::string const& file_name, po::variables_map& vm);

    static void luaopen(lua_State* L);

private:
    static void parse_command_line(po::command_line_parser& parser, po::variables_map& vm);

    /** module-independent options */
    po::options_description desc_;
    /** module-specific options */
    std::map<std::string, po::options_description> desc_module_;
    /** module namespaces */
    std::vector<std::string> sections_;
};

} // namespace halmd

#endif /* ! HALMD_OPTIONS_PARSER_HPP */
