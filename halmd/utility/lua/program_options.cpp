/*
 * Copyright © 2008-2012  Peter Colberg
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

#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/cmdline.hpp>
#include <luabind/luabind.hpp>
#include <luabind/adopt_policy.hpp>
#include <luabind/operator.hpp> // luabind::tostring
#include <luabind/return_reference_to_policy.hpp>
#include <luabind/shared_ptr_converter.hpp>
#include <stdint.h>

#include <halmd/config.hpp>
#include <halmd/utility/lua/vector_converter.hpp>

namespace po = boost::program_options;

using namespace boost;
using namespace std;

namespace std {

template <typename T>
static ostream& operator<<(ostream& os, vector<T> const& value)
{
    typename vector<T>::const_iterator i = value.begin();
    if (i != value.end()) {
        os << *i;
        for (++i; i != value.end(); ++i) {
            os << " " << *i;
        }
    }
    return os;
}

} // namespace std

namespace halmd {

template <typename T>
static po::typed_value<T>*
default_value(po::typed_value<T>* semantic, T const& value)
{
    return semantic->default_value(value);
}

template <typename T>
static po::typed_value<T>*
default_value_textual(po::typed_value<T>* semantic, T const& value, string const& textual)
{
    return semantic->default_value(value, textual);
}

template <typename T>
static po::typed_value<T>*
implicit_value(po::typed_value<T>* semantic, T const& value)
{
    return semantic->implicit_value(value);
}

template <typename T>
static po::typed_value<T>*
implicit_value_textual(po::typed_value<T>* semantic, T const& value, string const& textual)
{
    return semantic->implicit_value(value, textual);
}

template <typename T>
static void notify(luabind::object const& functor, T const& value)
{
    try {
        luabind::call_function<void>(functor, value);
    }
    catch (luabind::error const& e) {
        throw runtime_error(lua_tostring(e.state(), -1));
    }
}

template <typename T>
static po::typed_value<T>*
notifier(po::typed_value<T>* semantic, luabind::object const& functor)
{
    return semantic->notifier(bind(&notify<T>, functor, _1));
}

template <typename T>
struct typed_value_wrapper : po::typed_value<T>, luabind::wrap_base
{
    typed_value_wrapper() : po::typed_value<T>(0) {}
};

template <>
struct typed_value_wrapper<bool> : po::typed_value<bool>, luabind::wrap_base
{
    typed_value_wrapper() : po::typed_value<bool>(0)
    {
        default_value(false);
        zero_tokens();
    }
};

struct untyped_value_wrapper : po::untyped_value, luabind::wrap_base
{
    untyped_value_wrapper() : po::untyped_value(true) {} // zero tokens
};

static void
add_option(po::options_description& self, shared_ptr<po::option_description> desc)
{
    self.add(desc);
}

static void
add_options(po::options_description& self, po::options_description const& desc)
{
    self.add(desc);
}

static po::command_line_parser
disallow_guessing(po::command_line_parser& parser)
{
    using namespace boost::program_options::command_line_style;
    return parser.style(default_style & ~allow_guessing);
}

static void
variables_map_store(po::variables_map& vm, po::parsed_options const& options)
{
    po::store(options, vm);
}

static void
variables_map_notify(po::variables_map& vm)
{
    po::notify(vm);
}

template <typename T>
static void typed_value(lua_State* L, char const* name)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("program_options")
        [
            namespace_("value")
            [
                class_<po::typed_value<T>, typed_value_wrapper<T>, po::value_semantic>(name)
                    .def(constructor<>())
                    .def("default_value", &default_value<T>, return_reference_to(_1))
                    .def("default_value", &default_value_textual<T>, return_reference_to(_1))
                    .def("implicit_value", &implicit_value<T>, return_reference_to(_1))
                    .def("implicit_value", &implicit_value_textual<T>, return_reference_to(_1))
                    .def("notifier", &notifier<T>, return_reference_to(_1))
                    .def("required", &po::typed_value<T>::required, return_reference_to(_1))

            ]

          , namespace_("multivalue")
            [
                class_<po::typed_value<vector<T> >, typed_value_wrapper<vector<T> >, po::value_semantic>(name)
                    .def(constructor<>())
                    .def("default_value", &default_value<vector <T> >, return_reference_to(_1))
                    .def("default_value", &default_value_textual<vector <T> >, return_reference_to(_1))
                    .def("implicit_value", &implicit_value<vector <T> >, return_reference_to(_1))
                    .def("implicit_value", &implicit_value_textual<vector <T> >, return_reference_to(_1))
                    .def("notifier", &notifier<vector<T> >, return_reference_to(_1))
                    .def("required", &po::typed_value<vector<T> >::required, return_reference_to(_1))
                    .def("composing", &po::typed_value<vector<T> >::composing, return_reference_to(_1))
                    .def("multitoken", &po::typed_value<vector<T> >::multitoken, return_reference_to(_1))
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_utility_lua_program_options(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("program_options")
        [
            class_<po::value_semantic>("value_semantic")
                .property("name", &po::value_semantic::name)
                .property("min_tokens", &po::value_semantic::min_tokens)
                .property("max_tokens", &po::value_semantic::max_tokens)
                .property("is_composing", &po::value_semantic::is_composing)
                .property("is_required", &po::value_semantic::is_required)

          , class_<po::untyped_value, untyped_value_wrapper, po::value_semantic>("untyped_value")
                .def(constructor<>())

          , class_<po::option_description>("option_description")
                .def(constructor<char const*, po::value_semantic*>(), adopt(_3))
                .def(constructor<char const*, po::value_semantic*, char const*>(), adopt(_3))
                .def("key", &po::option_description::key)
                .property("description", &po::option_description::description)
                .property("long_name", &po::option_description::long_name)
                .property("semantic", &po::option_description::semantic)
                .property("format_name", &po::option_description::format_name)
                .property("format_parameter", &po::option_description::format_parameter)

          , class_<po::options_description>("options_description")
                .def(constructor<>())
                .def(constructor<unsigned int>())
                .def(constructor<unsigned int, unsigned int>())
                .def(constructor<string>())
                .def(constructor<string, unsigned int>())
                .def(constructor<string, unsigned int, unsigned int>())
                .def(tostring(const_self))
                .def("add", &add_option)
                .def("add", &add_options)
                .property("options", &po::options_description::options)

          , class_<po::positional_options_description>("positional_options_description")
                .def(constructor<>())
                .def("add", &po::positional_options_description::add)

          , class_<po::command_line_parser>("command_line_parser")
                .def(constructor<vector<string> const&>())
                .def("options", &po::command_line_parser::options, return_reference_to(_1))
                .def("positional", &po::command_line_parser::positional, return_reference_to(_1))
                .def("allow_unregistered", &po::command_line_parser::allow_unregistered, return_reference_to(_1))
                .def("disallow_guessing", &disallow_guessing, return_reference_to(_1))
                .def("run", &po::command_line_parser::run)

          , class_<po::parsed_options>("parsed_options")

          , class_<po::variables_map>("variables_map")
                .def(constructor<>())
                .def("store", &variables_map_store)
                .def("notify", &variables_map_notify)
                .def("count", &po::variables_map::count)
        ]
    ];
    typed_value<bool>(L, "boolean");
    typed_value<string>(L, "string");
    typed_value<int32_t>(L, "int32");
    typed_value<int64_t>(L, "int64");
    typed_value<uint32_t>(L, "uint32");
    typed_value<uint64_t>(L, "uint64");
    typed_value<float>(L, "float32");
    typed_value<double>(L, "float64");
    return 0;
}

} // namespace halmd
