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

#include <h5xx/h5xx.hpp>
#include <luabind/luabind.hpp>
#include <luabind/exception_handler.hpp>

#include <halmd/utility/lua_wrapper/hdf5.hpp>
#include <halmd/utility/lua_wrapper/long_long_converter.hpp> // *int64_t on x86
#include <halmd/utility/lua_wrapper/vector_converter.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace lua_wrapper { namespace hdf5
{

/**
 * Wrap C++ type for Luabind function overloading
 */
template <typename T>
struct type_wrapper {};

/**
 * Read attribute in HDF5 group/dataset if exists, otherwise return nil
 */
template <typename T>
static luabind::object read_attribute(
    H5::H5Object const& object, string const& name, type_wrapper<T> const&
  , lua_State* L
)
{
    if (h5xx::exists_attribute(object, name)) {
        T value = h5xx::read_attribute<T>(object, name);
        return luabind::object(L, cref(value));
    }
    return luabind::object(); // == nil
}

/**
 * Write attribute to HDF5 group/dataset
 */
template <typename T>
static void write_attribute(
    H5::H5Object const& object, string const& name, type_wrapper<T> const&
  , T const& value
)
{
    h5xx::write_attribute<T>(object, name, value);
}

/**
 * Translate HDF5 C++ exception to Lua error message
 *
 * H5::Exception does not derive from std::exception.
 */
static int translate_h5_exception(lua_State* L, H5::Exception const& e)
{
    lua_pushliteral(L, "HDF5 exception");
    return 1;
}

/**
 * Register HDF5 classes and functions with Lua
 */
int luaopen(lua_State* L)
{
    using namespace luabind;
    register_exception_handler<H5::Exception>(&translate_h5_exception);
    module(L, "libhalmd")
    [
        namespace_("h5")
        [
            class_<H5::IdComponent>("id")

          , class_<H5::AbstractDs>("abstract_dataset")

          , class_<H5::CommonFG>("common_fg")
                .def("open_group", &h5xx::open_group)
                .def("exists_group", &h5xx::exists_group)

          , class_<H5::H5Object, H5::IdComponent>("object")
                .def("read_attribute", &read_attribute<bool>)
                .def("read_attribute", &read_attribute<int>)
                .def("read_attribute", &read_attribute<unsigned int>)
                .def("read_attribute", &read_attribute<int64_t>)
                .def("read_attribute", &read_attribute<uint64_t>)
                .def("read_attribute", &read_attribute<double>)
                .def("read_attribute", &read_attribute<string>)
                .def("read_attribute", &read_attribute<vector<int> >)
                .def("read_attribute", &read_attribute<vector<unsigned int> >)
                .def("read_attribute", &read_attribute<vector<int64_t> >)
                .def("read_attribute", &read_attribute<vector<uint64_t> >)
                .def("read_attribute", &read_attribute<vector<double> >)
                .def("write_attribute", &write_attribute<bool>)
                .def("write_attribute", &write_attribute<int>)
                .def("write_attribute", &write_attribute<unsigned int>)
                .def("write_attribute", &write_attribute<int64_t>)
                .def("write_attribute", &write_attribute<uint64_t>)
                .def("write_attribute", &write_attribute<double>)
                .def("write_attribute", &write_attribute<string>)
                .def("write_attribute", &write_attribute<vector<int> >)
                .def("write_attribute", &write_attribute<vector<unsigned int> >)
                .def("write_attribute", &write_attribute<vector<int64_t> >)
                .def("write_attribute", &write_attribute<vector<uint64_t> >)
                .def("write_attribute", &write_attribute<vector<double> >)

          , class_<H5::H5File, bases<H5::IdComponent, H5::CommonFG> >("file")

          , class_<H5::Group, bases<H5::H5Object, H5::CommonFG> >("group")

          , class_<H5::DataSet, bases<H5::H5Object, H5::AbstractDs> >("dataset")
                .def("exists_dataset", &h5xx::exists_dataset)

          , class_<type_wrapper<bool> >("bool")
                .def(constructor<>())
          , class_<type_wrapper<int> >("int")
                .def(constructor<>())
          , class_<type_wrapper<unsigned int> >("uint")
                .def(constructor<>())
          , class_<type_wrapper<int64_t> >("int64")
                .def(constructor<>())
          , class_<type_wrapper<uint64_t> >("uint64")
                .def(constructor<>())
          , class_<type_wrapper<double> >("float")
                .def(constructor<>())
          , class_<type_wrapper<string> >("string")
                .def(constructor<>())

          , class_<type_wrapper<vector<int> > >("int_array")
                .def(constructor<>())
          , class_<type_wrapper<vector<unsigned int> > >("uint_array")
                .def(constructor<>())
          , class_<type_wrapper<vector<int64_t> > >("int64_array")
                .def(constructor<>())
          , class_<type_wrapper<vector<uint64_t> > >("uint64_array")
                .def(constructor<>())
          , class_<type_wrapper<vector<double> > >("float_array")
                .def(constructor<>())
        ]
    ];
    return 0;
}

}} // namespace lua_wrapper::hdf5

} // namespace halmd
