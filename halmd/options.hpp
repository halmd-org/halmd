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

#ifndef HALMD_OPTIONS_HPP
#define HALMD_OPTIONS_HPP

#include <boost/array.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/program_options.hpp>
#include <lua.hpp>
#include <string>
#include <vector>

#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

namespace halmd
{

// abbreviated namespace with top-level scope
namespace po = boost::program_options;

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

/**
 * Accumulating program option value
 */
template <typename T, typename charT = char>
class accumulating_value
  : public boost::program_options::typed_value<T,charT>
{
    //
    // Originally written by Bryan Green <bgreen@nas.nasa.gov>
    // http://article.gmane.org/gmane.comp.lib.boost.user/29084
    //

public:
    accumulating_value(T* store_to = 0)
      : boost::program_options::typed_value<T,charT>(store_to)
      , origin(0)
    {
        (void) boost::program_options::typed_value<T,charT>::zero_tokens();
    }

    accumulating_value* default_value(const T &v)
    {
        // setting a default value sets the origin to that value
        origin = v;
        (void) boost::program_options::typed_value<T,charT>::default_value(v);
        return this;
    }

    accumulating_value* default_value(const T &v, const std::string& textual)
    {
        // setting a default value sets the origin to that value
        origin = v;
        (void) boost::program_options::typed_value<T,charT>::default_value(v, textual);
        return this;
    }

    void xparse(boost::any& value_store, const std::vector<std::basic_string<charT> >& new_tokens) const
    {
        // if this is the first occurrence of the option, initialize it to the origin
        if (value_store.empty()) {
            value_store = boost::any(origin);
        }
        ++boost::any_cast<T&>(value_store);
    }

private:
    T origin; //< the numeric origin from which to increment upward
};

template <typename T>
accumulating_value<T>* accum_value(T *v)
{
    accumulating_value<T>* r = new accumulating_value<T>(v);
    return r;
}

template <typename T>
accumulating_value<T>* accum_value()
{
    return accum_value<T>(0);
}

} // namespace halmd

namespace std
{

/**
 * extract comma-separated option values into fixed-size array
 */
template <typename T, size_t size>
std::istream& operator>>(std::istream& is, boost::array<T, size>& value)
{
    BOOST_FOREACH(T& v, value) {
        std::string str;
        getline(is, str, ',');
        v = boost::lexical_cast<T>(str);
    }
    return is;
}

template <typename T, size_t size>
std::ostream& operator<<(std::ostream& os, boost::array<T, size> const& value)
{
    BOOST_FOREACH(T const& v, value) {
        if (&v != &value.front()) {
            os << ',';
        }
        os << v;
    }
    return os;
}

/**
 * extract comma-separated option values into variable-size array
 */
template <typename T>
std::istream& operator>>(std::istream& is, boost::multi_array<T, 1>& value)
{
    std::vector<T> v;
    std::string str;
    while (!is.eof()) {
        getline(is, str, ',');
        v.push_back(boost::lexical_cast<T>(str));
    }
    boost::array<size_t, 1> extents = boost::assign::list_of(v.size());
    value.resize(extents);
    std::copy(v.begin(), v.end(), value.begin());
    return is;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, boost::multi_array<T, 1> const& value)
{
    BOOST_FOREACH(T const& v, value) {
        if (&v != &(*value.begin())) {
            os << ',';
        }
        os << v;
    }
    return os;
}

} // namespace std

namespace luabind
{

/**
 * Luabind converter for Boost.Program_options variables_map
 */
template <>
struct default_converter<boost::program_options::variables_map>
  : native_converter_base<boost::program_options::variables_map>
{
    //! convert from C++ to Lua
    void to(lua_State* L, boost::program_options::variables_map const& vm)
    {
        luabind::object table = luabind::newtable(L);
        boost::program_options::variables_map::const_iterator it, end = vm.end();
        for (it = vm.begin(); it != end; ++it) {
            table[it->first] = it->second;
        }
        table.push(L);
    }
};

template <>
struct default_converter<boost::program_options::variables_map const&>
  : default_converter<boost::program_options::variables_map> {};

} // namespace luabind

#endif /* ! HALMD_OPTIONS_HPP */
