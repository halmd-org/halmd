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

#include <halmd/io/logger.hpp> //< logger::warning
#include <halmd/utility/options_parser.hpp>
#include <halmd/utility/date_time.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
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
        ("verbose,v", po::accum_value<int>()->default_value(logger::warning),
         "increase verbosity")
        ("version",
         "output version and exit")
        ("help",
         "display this help and exit")
        ;
}

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
    // disallow abbreviated options, which breaks forward compatibility of
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
 * print options parser error message to stderr
 *
 * @param error exception deriving from std::exception
 */
void options_parser::print_error(std::exception const& error) const
{
    cerr << PROGRAM_NAME ": " << error.what() << endl;
    cerr << "Try `" PROGRAM_NAME " --help' for more information." << endl;
}

/**
 * print options help message to stdout
 */
void options_parser::print_help() const
{
    cout << "Usage: " PROGRAM_NAME " [OPTION]..." << endl << endl
         << desc_ << endl;
}

/**
 * print version information to stdout
 */
void options_parser::print_version() const
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

template <typename T>
static po::extended_typed_value<T>* po_value()
{
    return po::value<T>();
}

static po::extended_typed_value<bool>* po_bool_switch()
{
    return po::bool_switch();
}

template <typename T>
static void po_call_notifier(luabind::object const& f, T const& value)
{
    luabind::call_function<void>(f, cref(value));
}

template <typename T>
static po::extended_typed_value<T>* po_notifier(
    po::extended_typed_value<T>* v, luabind::object const& f
)
{
    return v->notifier(bind(&po_call_notifier<T>, f, _1));
}

static void po_add_option_description(
    po::options_description& desc, char const* name
  , po::value_semantic const* semantic, char const* description
)
{
    desc.add_options()(name, semantic, description);
}

static void po_add_options_description(
    po::options_description& desc, po::options_description const& other
)
{
    desc.add(other);
}

/**
 * register Lua C++ wrapper
 */
void options_parser::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "halmd_wrapper")
    [
        namespace_("po")
        [
            class_<po::options_description>("options_description")
                .def(constructor<>())
                .def(constructor<string>())
                .def("add", &po_add_option_description)
                .def("add", &po_add_options_description)

          , class_<po::variable_value>("variable_value")
                .def(constructor<>())
                .def("empty", &po::variable_value::empty)
                .def("defaulted", &po::variable_value::defaulted)
                .def("value", (any const& (po::variable_value::*)() const) &po::variable_value::value)
                //< only return-by-value is supported by Luabind boost::any converter

          , class_<po::value_semantic>("value_semantic")

          , class_<po::extended_typed_value<bool>, po::value_semantic>("typed_value_bool")
                .def("notifier", &po_notifier<bool>)
                .def("conflicts", &po::extended_typed_value<bool>::conflicts)
                .def("depends", &po::extended_typed_value<bool>::depends)

          , class_<po::extended_typed_value<int>, po::value_semantic>("typed_value_int")
                .def("notifier", &po_notifier<int>)
                .def("conflicts", &po::extended_typed_value<int>::conflicts)
                .def("depends", &po::extended_typed_value<int>::depends)

          , class_<po::extended_typed_value<unsigned int>, po::value_semantic>("typed_value_uint")
                .def("notifier", &po_notifier<unsigned int>)
                .def("conflicts", &po::extended_typed_value<unsigned int>::conflicts)
                .def("depends", &po::extended_typed_value<unsigned int>::depends)

          , class_<po::extended_typed_value<int64_t>, po::value_semantic>("typed_value_int64")
                .def("notifier", &po_notifier<int64_t>)
                .def("conflicts", &po::extended_typed_value<int64_t>::conflicts)
                .def("depends", &po::extended_typed_value<int64_t>::depends)

          , class_<po::extended_typed_value<uint64_t>, po::value_semantic>("typed_value_uint64")
                .def("notifier", &po_notifier<uint64_t>)
                .def("conflicts", &po::extended_typed_value<uint64_t>::conflicts)
                .def("depends", &po::extended_typed_value<uint64_t>::depends)

          , class_<po::extended_typed_value<double>, po::value_semantic>("typed_value_float")
                .def("notifier", &po_notifier<double>)
                .def("conflicts", &po::extended_typed_value<double>::conflicts)
                .def("depends", &po::extended_typed_value<double>::depends)

          , class_<po::extended_typed_value<string>, po::value_semantic>("typed_value_string")
                .def("notifier", &po_notifier<string>)
                .def("conflicts", &po::extended_typed_value<string>::conflicts)
                .def("depends", &po::extended_typed_value<string>::depends)

          , def("bool_switch", &po_bool_switch)
          , def("int", &po_value<int>)
          , def("uint", &po_value<unsigned int>)
          , def("int64", &po_value<int64_t>)
          , def("uint64", &po_value<uint64_t>)
          , def("float", &po_value<double>)
          , def("string", &po_value<string>)

          , class_<po::extended_typed_value<multi_array<int, 1> >, po::value_semantic>("typed_value_int_array")
                .def("notifier", &po_notifier<int>)
                .def("conflicts", &po::extended_typed_value<multi_array<int, 1> >::conflicts)
                .def("depends", &po::extended_typed_value<multi_array<int, 1> >::depends)

          , class_<po::extended_typed_value<multi_array<unsigned int, 1> >, po::value_semantic>("typed_value_uint_array")
                .def("notifier", &po_notifier<unsigned int>)
                .def("conflicts", &po::extended_typed_value<multi_array<unsigned int, 1> >::conflicts)
                .def("depends", &po::extended_typed_value<multi_array<unsigned int, 1> >::depends)

          , class_<po::extended_typed_value<multi_array<int64_t, 1> >, po::value_semantic>("typed_value_int64_array")
                .def("notifier", &po_notifier<int64_t>)
                .def("conflicts", &po::extended_typed_value<multi_array<int64_t, 1> >::conflicts)
                .def("depends", &po::extended_typed_value<multi_array<int64_t, 1> >::depends)

          , class_<po::extended_typed_value<multi_array<uint64_t, 1> >, po::value_semantic>("typed_value_uint64_array")
                .def("notifier", &po_notifier<uint64_t>)
                .def("conflicts", &po::extended_typed_value<multi_array<uint64_t, 1> >::conflicts)
                .def("depends", &po::extended_typed_value<multi_array<uint64_t, 1> >::depends)

          , class_<po::extended_typed_value<multi_array<double, 1> >, po::value_semantic>("typed_value_float_array")
                .def("notifier", &po_notifier<double>)
                .def("conflicts", &po::extended_typed_value<multi_array<double, 1> >::conflicts)
                .def("depends", &po::extended_typed_value<multi_array<double, 1> >::depends)

          , def("int_array", &po_value<multi_array<int, 1> >)
          , def("uint_array", &po_value<multi_array<unsigned int, 1> >)
          , def("int64_array", &po_value<multi_array<int64_t, 1> >)
          , def("uint64_array", &po_value<multi_array<uint64_t, 1> >)
          , def("float_array", &po_value<multi_array<double, 1> >)
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &options_parser::luaopen
    ];
}

} // namespace halmd
