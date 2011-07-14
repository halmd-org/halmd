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

#include <iomanip>
#include <luabind/luabind.hpp>
#include <sstream>
#include <stdint.h> // <cstdint> is C++0x

#include <halmd/config.hpp>
#include <halmd/utility/lua/array_converter.hpp>
#include <halmd/utility/lua/long_long_converter.hpp> // *int64_t on x86
#include <halmd/utility/lua/map_converter.hpp>
#include <halmd/utility/program_options/program_options.hpp>

using namespace boost;
using namespace std;

namespace halmd {

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
static void po_lua_notifier(luabind::object const& notifier, T const& value)
{
    using namespace luabind;

    object result;
    try {
        result = call_function<object>(notifier, value);
    }
    catch (error const& e) {
        throw po::error(lua_tostring(e.state(), -1));
    }

    if (result) {
        const_cast<T&>(value) = object_cast<T>(result);
    }
}

template <typename T>
static po::extended_typed_value<T>* po_notifier(
    po::extended_typed_value<T>* v, luabind::object const& notifier
)
{
    return v->notifier(bind(&po_lua_notifier<T>, notifier, _1));
}

template <typename T>
static void po_choices_notifier(map<T, string> const& choices, T const& value)
{
    if (choices.find(value) == choices.end()) {
        typename map<T, string>::const_iterator i, ie;
        size_t pad = 21; //< minimum padding, equivalent to --help output
        for (i = choices.begin(), ie = choices.end(); i != ie; ++i) {
            pad = max(pad, i->first.size());
        }

        stringstream s;
        s << "invalid option value '" << value << "'" << endl << endl;

        s << "The choices for the option are:" << endl;
        for (i = choices.begin(), ie = choices.end(); i != ie; ++i) {
            s << "  " << left << setw(pad) << i->first << " " << i->second << endl;
        }

        throw po::error(s.str());
    }
}

template <typename T>
static po::extended_typed_value<T>* po_choices(
    po::extended_typed_value<T>* v, map<T, string> const& choices
)
{
    return v->notifier(bind(&po_choices_notifier<T>, choices, _1));
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
 * register Boost Program_otions with Lua
 */
HALMD_LUA_API int luaopen_libhalmd_utility_lua_program_options(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("po")
        [
            class_<po::value_semantic>("value_semantic")

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
                .def("choices", &po_choices<string>)

          , class_<po::extended_typed_value<multi_array<int, 1> >, po::value_semantic>("typed_value_int_array")
                .def("notifier", &po_notifier<multi_array<int, 1> >)
                .def("conflicts", &po::extended_typed_value<multi_array<int, 1> >::conflicts)
                .def("depends", &po::extended_typed_value<multi_array<int, 1> >::depends)

          , class_<po::extended_typed_value<multi_array<unsigned int, 1> >, po::value_semantic>("typed_value_uint_array")
                .def("notifier", &po_notifier<multi_array<unsigned int, 1> >)
                .def("conflicts", &po::extended_typed_value<multi_array<unsigned int, 1> >::conflicts)
                .def("depends", &po::extended_typed_value<multi_array<unsigned int, 1> >::depends)

          , class_<po::extended_typed_value<multi_array<int64_t, 1> >, po::value_semantic>("typed_value_int64_array")
                .def("notifier", &po_notifier<multi_array<int64_t, 1> >)
                .def("conflicts", &po::extended_typed_value<multi_array<int64_t, 1> >::conflicts)
                .def("depends", &po::extended_typed_value<multi_array<int64_t, 1> >::depends)

          , class_<po::extended_typed_value<multi_array<uint64_t, 1> >, po::value_semantic>("typed_value_uint64_array")
                .def("notifier", &po_notifier<multi_array<uint64_t, 1> >)
                .def("conflicts", &po::extended_typed_value<multi_array<uint64_t, 1> >::conflicts)
                .def("depends", &po::extended_typed_value<multi_array<uint64_t, 1> >::depends)

          , class_<po::extended_typed_value<multi_array<double, 1> >, po::value_semantic>("typed_value_float_array")
                .def("notifier", &po_notifier<multi_array<double, 1> >)
                .def("conflicts", &po::extended_typed_value<multi_array<double, 1> >::conflicts)
                .def("depends", &po::extended_typed_value<multi_array<double, 1> >::depends)

          , def("bool_switch", &po_bool_switch)
          , def("int", &po_value<int>)
          , def("uint", &po_value<unsigned int>)
          , def("int64", &po_value<int64_t>)
          , def("uint64", &po_value<uint64_t>)
          , def("float", &po_value<double>)
          , def("string", &po_value<string>)
          , def("int_array", &po_value<multi_array<int, 1> >)
          , def("uint_array", &po_value<multi_array<unsigned int, 1> >)
          , def("int64_array", &po_value<multi_array<int64_t, 1> >)
          , def("uint64_array", &po_value<multi_array<uint64_t, 1> >)
          , def("float_array", &po_value<multi_array<double, 1> >)

          , class_<po::options_description>("options_description")
                .def(constructor<>())
                .def(constructor<string>())
                .def("add", &po_add_option_description)
                .def("add", &po_add_options_description)
        ]
    ];
    return 0;
}

} // namespace halmd
