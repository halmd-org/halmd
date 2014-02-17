/*
 * Copyright © 2014      Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#include <algorithm>
#include <boost/algorithm/string/find_iterator.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/any.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/range/iterator_range.hpp> // boost::copy_range
#include <boost/utility/enable_if.hpp>
#include <luaponte/luaponte.hpp>
#include <luaponte/adopt_policy.hpp>
#include <luaponte/operator.hpp> // luaponte::tostring
#include <luaponte/return_reference_to_policy.hpp>
#include <luaponte/shared_ptr_converter.hpp>
#include <sstream>
#include <stdint.h>

#include <halmd/config.hpp>
#include <halmd/numeric/cast.hpp>
#include <halmd/utility/lua/ublas.hpp>
#include <halmd/utility/lua/vector_converter.hpp>
#include <halmd/utility/program_options.hpp>

namespace po = boost::program_options;
namespace ublas = boost::numeric::ublas;

using namespace boost::algorithm; // first_finder, make_split_iterator, split_iterator
using namespace std;

namespace std {

/**
 * Ensure arithmetic option value is exactly representable as Lua number.
 *
 * Note that floating-point promotion from float to double (which is the
 * default for lua_Number) is exact, refer to ISO/IEC 14882:1998, 4.6.1.
 * “An rvalue of type float can be converted to an rvalue of type double.
 * The value is unchanged.”
 */
template <typename T>
static typename boost::enable_if<boost::is_arithmetic<T>, void>::type
validate(boost::any& v, vector<string> const& values, T*, int)
{
    po::validators::check_first_occurrence(v);
    string s(po::validators::get_single_string(values));
    try {
        T value = boost::lexical_cast<T>(s);
        halmd::checked_narrowing_cast<lua_Number>(value);
        v = boost::any(value);
    }
    catch (exception const&) {
        throw po::invalid_option_value(s);
    }
}

/**
 * Read Boost uBLAS vector from input stream
 *
 * A vector is represented as a comma-delimited string, e.g. 1,2,3,4
 */
template <typename T>
void validate(boost::any& v, vector<string> const& values, ublas::vector<T>*, int)
{
    po::validators::check_first_occurrence(v);
    string s(po::validators::get_single_string(values));
    ublas::vector<T> value;
    vector<T> element;
    for (split_iterator<string::iterator> i = make_split_iterator(s, first_finder(",", is_equal()));
         i != split_iterator<string::iterator>();
         ++i)
    {
        boost::any v;
        vector<string> values;
        values.push_back(boost::copy_range<string>(*i));
        validate(v, values, (T*)0, 0);
        element.push_back(boost::any_cast<T>(v));
    }
    value.resize(element.size());
    copy(element.begin(), element.end(), value.begin());
    v = value;
}

/**
 * Read Boost uBLAS matrix from input stream
 *
 * A matrix is represented as colon-delimited rows with each row as a
 * comma-delimited string, e.g. 11,12,13:21,22,23:31,32,33
 */
template <typename T>
void validate(boost::any& v, vector<string> const& values, ublas::matrix<T>*, int)
{
    po::validators::check_first_occurrence(v);
    string s(po::validators::get_single_string(values));
    ublas::matrix<T> value;
    vector<ublas::vector<T> > row;
    for (split_iterator<string::iterator> i = make_split_iterator(s, first_finder(":", is_equal()));
         i != split_iterator<string::iterator>();
         ++i)
    {
        boost::any v;
        vector<string> values;
        values.push_back(boost::copy_range<string>(*i));
        validate(v, values, (ublas::vector<T>*)0, 0);
        row.push_back(boost::any_cast<ublas::vector<T> >(v));
    }
    if (!row.empty()) {
        if (row.front().size() == 1) {
            ublas::symmetric_matrix<T, ublas::lower> m(row.size(), row.size());
            for (size_t i = 0; i < row.size(); ++i) {
                if (!(row[i].size() == i + 1)) {
                    throw po::invalid_option_value(s);
                }
                copy(row[i].begin(), row[i].end(), ublas::row(m, i).begin());
            }
            value = m;
        }
        else if (row.back().size() == 1) {
            ublas::symmetric_matrix<T, ublas::upper> m(row.size(), row.size());
            for (size_t i = 0; i < row.size(); ++i) {
                if (!(row[i].size() == (m.size2() - i))) {
                    throw po::invalid_option_value(s);
                }
                copy(row[i].begin(), row[i].end(), ublas::row(m, i).begin());
            }
            value = m;
        }
        else {
            value.resize(row.size(), row.front().size());
            for (size_t i = 0; i < row.size(); ++i) {
                if (!(row[i].size() == value.size2())) {
                    throw po::invalid_option_value(s);
                }
                ublas::row(value, i) = row[i];
            }
        }
    }
    v = value;
}

/**
 * Write Boost uBLAS vector to output stream
 */
template <typename T>
ostream& operator<<(ostream& os, ublas::vector<T> const& value)
{
    for (size_t i = 0; i < value.size(); ++i) {
        if (i > 0) {
            os << ',';
        }
        os << value(i);
    }
    return os;
}

/**
 * Write Boost uBLAS matrix to output stream
 */
template <typename T>
ostream& operator<<(ostream& os, ublas::matrix<T> const& value)
{
    for (size_t i = 0; i < value.size1(); ++i) {
        if (i > 0) {
            os << ':';
        }
        ublas::vector<T> const& row = ublas::matrix_row<ublas::matrix<T> const>(value, i);
        os << row;
    }
    return os;
}

/**
 * Write STL vector to output stream
 */
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

// boost::program_options uses boost::lexical_cast to format default values for
// output, which messes up with rounding of floating-point numbers, e.g., 0.1
// becomes "0.10000000000000001". To resolve the issue, we specialise the
// template for T=float,double and let std::ostringstream convert the value.
template <>
po::typed_value<float>*
default_value(po::typed_value<float>* semantic, float const& value)
{
    std::ostringstream ss;
    ss << value;
    return semantic->default_value(value, ss.str());
}

template <>
po::typed_value<double>*
default_value(po::typed_value<double>* semantic, double const& value)
{
    std::ostringstream ss;
    ss << value;
    return semantic->default_value(value, ss.str());
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
static void notify(luaponte::object const& functor, T const& value)
{
    luaponte::object args(luaponte::from_stack(functor.interpreter(), -1));
    try {
        luaponte::call_function<void>(functor, args, value);
    }
    catch (luaponte::error const& e) {
        string error(lua_tostring(e.state(), -1));
        lua_pop(e.state(), 1);
        throw runtime_error(error);
    }
}

template <typename T>
static po::typed_value<T>*
notifier(po::typed_value<T>* semantic, luaponte::object const& functor)
{
    return semantic->notifier([=](T const& value) {
        notify(functor, value);
    });
}

template <typename T>
struct typed_value_wrapper : po::typed_value<T>, luaponte::wrap_base
{
    typed_value_wrapper() : po::typed_value<T>(0) {}
};

template <>
struct typed_value_wrapper<bool> : po::typed_value<bool>, luaponte::wrap_base
{
    typed_value_wrapper() : po::typed_value<bool>(0)
    {
        default_value(false);
        zero_tokens();
    }
};

struct untyped_value_wrapper : po::untyped_value, luaponte::wrap_base
{
    untyped_value_wrapper() : po::untyped_value(true) {} // zero tokens
};

template <typename T>
static accumulating_value<T>*
accum_default_value(accumulating_value<T>* semantic, T const& value)
{
    return semantic->default_value(value);
}

template <typename T>
static accumulating_value<T>*
accum_default_value_textual(accumulating_value<T>* semantic, T const& value, string const& textual)
{
    return semantic->default_value(value, textual);
}

template <typename T>
struct accum_value_wrapper : accumulating_value<T>, luaponte::wrap_base
{
    accum_value_wrapper() : accumulating_value<T>(0) {}
};

static void
add_option(po::options_description& self, boost::shared_ptr<po::option_description> desc)
{
    self.add(desc);
}

static void
add_options(po::options_description& self, po::options_description const& desc)
{
    self.add(desc);
}

static po::command_line_parser&
disallow_guessing(po::command_line_parser& parser)
{
    using namespace boost::program_options::command_line_style;
    return parser.style(default_style & ~allow_guessing);
}

static pair<string, string>
split_argument(string const& arg)
{
    size_t pos = arg.find('=');
    if (pos != string::npos) {
        return make_pair(arg.substr(0, pos), arg.substr(pos + 1));
    }
    return make_pair(arg, string());
}

static po::command_line_parser&
group_parser(po::command_line_parser& parser)
{
    return parser.extra_parser(&split_argument);
}

static void
variables_map_store(po::variables_map& vm, po::parsed_options const& options)
{
    po::store(options, vm);
}

struct scoped_push
{
    lua_State* const L;

    scoped_push(luaponte::object const& object) : L(object.interpreter())
    {
        object.push(L);
    }

    ~scoped_push()
    {
        lua_pop(L, 1);
    }
};

static luaponte::object
variables_map_notify(lua_State* L, po::variables_map& vm, luaponte::object const& args)
{
    scoped_push p(args);
    po::notify(vm);
    return args;
}

template <typename T>
static void typed_value(lua_State* L, char const* name, char const* value, char const* multi_value)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("program_options")
        [
            namespace_(value)
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

          , namespace_(multi_value)
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
    using namespace luaponte;
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

          , class_<accumulating_value<int>, accum_value_wrapper<int>, po::value_semantic>("accum_value")
                .def(constructor<>())
                .def("default_value", &accum_default_value<int>, return_reference_to(_1))
                .def("default_value", &accum_default_value_textual<int>, return_reference_to(_1))
                .def("notifier", &notifier<int>, return_reference_to(_1))
                .def("required", &po::typed_value<int>::required, return_reference_to(_1))

          , class_<po::option_description, boost::shared_ptr<po::option_description> >("option_description")
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
                .def("group_parser", &group_parser, return_reference_to(_1))
                .def("run", &po::command_line_parser::run)

          , class_<po::parsed_options>("parsed_options")

          , class_<po::variables_map>("variables_map")
                .def(constructor<>())
                .def("store", &variables_map_store)
                .def("notify", &variables_map_notify)
                .def("count", &po::variables_map::count)
        ]
    ];

    typed_value<bool                    >(L, "boolean", "value"       , "multi_value"       );
    typed_value<string                  >(L, "string" , "value"       , "multi_value"       );
    typed_value<int32_t                 >(L, "int32"  , "value"       , "multi_value"       );
    typed_value<int64_t                 >(L, "int64"  , "value"       , "multi_value"       );
    typed_value<uint32_t                >(L, "uint32" , "value"       , "multi_value"       );
    typed_value<uint64_t                >(L, "uint64" , "value"       , "multi_value"       );
    typed_value<float                   >(L, "float32", "value"       , "multi_value"       );
    typed_value<double                  >(L, "float64", "value"       , "multi_value"       );

    typed_value<ublas::vector<int32_t>  >(L, "int32"  , "value_vector", "multi_value_vector");
    typed_value<ublas::vector<int64_t>  >(L, "int64"  , "value_vector", "multi_value_vector");
    typed_value<ublas::vector<uint32_t> >(L, "uint32" , "value_vector", "multi_value_vector");
    typed_value<ublas::vector<uint64_t> >(L, "uint64" , "value_vector", "multi_value_vector");
    typed_value<ublas::vector<float>    >(L, "float32", "value_vector", "multi_value_vector");
    typed_value<ublas::vector<double>   >(L, "float64", "value_vector", "multi_value_vector");

    typed_value<ublas::matrix<int32_t>  >(L, "int32"  , "value_matrix", "multi_value_matrix");
    typed_value<ublas::matrix<int64_t>  >(L, "int64"  , "value_matrix", "multi_value_matrix");
    typed_value<ublas::matrix<uint32_t> >(L, "uint32" , "value_matrix", "multi_value_matrix");
    typed_value<ublas::matrix<uint64_t> >(L, "uint64" , "value_matrix", "multi_value_matrix");
    typed_value<ublas::matrix<float>    >(L, "float32", "value_matrix", "multi_value_matrix");
    typed_value<ublas::matrix<double>   >(L, "float64", "value_matrix", "multi_value_matrix");

    return 0;
}

} // namespace halmd
