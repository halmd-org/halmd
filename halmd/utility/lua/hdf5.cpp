/*
 * Copyright © 2014 Felix Höfling
 * Copyright © 2010 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <h5xx/h5xx.hpp>
#include <luaponte/luaponte.hpp>
#include <luaponte/exception_handler.hpp>
#include <luaponte/out_value_policy.hpp>

#include <halmd/config.hpp>
#include <halmd/utility/lua/vector_converter.hpp>

using namespace std;

namespace halmd {

/**
 * Wrap C++ type for Luabind function overloading
 */
template <typename T>
struct type_wrapper {};

/**
 * Read attribute in HDF5 group/dataset if exists, otherwise return nil
 */
template <typename T>
static luaponte::object read_attribute(
    H5::H5Object const& object, string const& name, type_wrapper<T> const&
  , lua_State* L
)
{
    if (h5xx::exists_attribute(object, name)) {
        T value = h5xx::read_attribute<T>(object, name);
        return luaponte::object(L, boost::cref(value));
    }
    return luaponte::object(); // == nil
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
 * Returns shape of dataset.
 */
static void wrap_shape(H5::DataSet const& dataset, std::vector<hsize_t>& shape)
{
    H5::DataSpace const& space = dataset.getSpace();
    shape.resize(space.getSimpleExtentNdims());
    space.getSimpleExtentDims(&*shape.begin());
}

/**
 * Returns a string encoding the data type of the dataset.
 */
static std::string wrap_type(H5::DataSet const& dataset)
{
    switch(dataset.getTypeClass())
    {
        case H5T_INTEGER:
            switch (dataset.getIntType().getSign())
            {
                case H5T_SGN_NONE:
                    return "unsigned_int";
                case H5T_SGN_2:
                    return "int";
                default:
                    throw std::runtime_error("unsupported dataset type: invalid integer sign");
            }
        case H5T_FLOAT:
#ifndef USE_HOST_SINGLE_PRECISION
            return "double";
#else
            return "float";
#endif
        case H5T_TIME:
            throw std::runtime_error("unsupported dataset type: time");
        case H5T_STRING:
            throw std::runtime_error("unsupported dataset type: string");
        case H5T_BITFIELD:
            throw std::runtime_error("unsupported dataset type: bitfield");
        case H5T_OPAQUE:
            throw std::runtime_error("unsupported dataset type: opaque");
        case H5T_COMPOUND:
            throw std::runtime_error("unsupported dataset type: compound");
        case H5T_REFERENCE:
            throw std::runtime_error("unsupported dataset type: reference");
        case H5T_ENUM:
            throw std::runtime_error("unsupported dataset type: enum");
        case H5T_VLEN:
            throw std::runtime_error("unsupported dataset type: vlen");
        case H5T_ARRAY:
            throw std::runtime_error("unsupported dataset type: array");
        default:
            throw std::runtime_error("unsupported dataset type: unknown");
    }
}

#if (H5_VERS_MAJOR > 1) || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 10)
#define HALMD_H5_COMMON_BASE H5::H5Location
/**
 * Open HDF5 DataSet.
 */
static H5::DataSet wrap_open_dataset(H5::H5Location const& group, std::string const& name)
{
    return group.openDataSet(name);
}

static bool wrap_exists_group(H5::H5Location const& loc, std::string const& name)
{
    using namespace h5xx;
    hid_t hid;
    H5E_BEGIN_TRY {
        hid = H5Gopen(loc.getId(), name.c_str(), H5P_DEFAULT);
        if (hid > 0) {
            H5Gclose(hid);
        }
    } H5E_END_TRY
    return (hid > 0);
}

static H5::Group wrap_open_group(H5::H5Location const& loc, std::string const& path)
{
    using namespace h5xx;
    hid_t group_id;
    H5E_BEGIN_TRY {
        group_id = H5Gopen(loc.getId(), path.c_str(), H5P_DEFAULT);
    } H5E_END_TRY
    if (group_id < 0) {
        H5::PropList pl = create_intermediate_group_property();
        group_id = H5Gcreate(loc.getId(), path.c_str(), pl.getId(), H5P_DEFAULT, H5P_DEFAULT);
    }
    if (group_id < 0) {
        throw error("failed to create group \"" + path + "\"");
    }
    return H5::Group(group_id);
}
#else
#define HALMD_H5_COMMON_BASE H5::CommonFG
/**
 * Open HDF5 DataSet.
 */
static H5::DataSet wrap_open_dataset(H5::CommonFG const& group, std::string const& name)
{
    return group.openDataSet(name);
}
#endif

/**
 * Register HDF5 classes and functions with Lua
 */
HALMD_LUA_API int luaopen_libhalmd_utility_lua_hdf5(lua_State* L)
{
    using namespace luaponte;
    register_exception_handler<H5::Exception>(&translate_h5_exception);
    module(L, "libhalmd")
    [
        namespace_("h5")
        [
            class_<H5::IdComponent>("id")

          , class_<H5::AbstractDs>("abstract_dataset")

#if (H5_VERS_MAJOR > 1) || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 10)
          , class_<H5::H5Location>("location")
                .def("open_dataset", &wrap_open_dataset)
                .def("open_group", &wrap_open_group)
                .def("exists_group", &wrap_exists_group)
#else
          , class_<H5::CommonFG>("common_fg")
                .def("open_dataset", &wrap_open_dataset)
                .def("open_group", &h5xx::open_group)
                .def("exists_group", &h5xx::exists_group)
#endif

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
                .def("read_attribute", &read_attribute<vector<string> >)

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
                .def("write_attribute", &write_attribute<vector<string> >)

          , class_<H5::H5File, bases<H5::IdComponent, HALMD_H5_COMMON_BASE> >("file")

          , class_<H5::Group, bases<H5::H5Object, HALMD_H5_COMMON_BASE> >("group")
                .def("exists_dataset", &h5xx::exists_dataset)

          , class_<H5::DataSet, bases<H5::H5Object, H5::AbstractDs> >("dataset")
                .property("shape", &wrap_shape, pure_out_value(_2))
                .property("type", &wrap_type)

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
          , class_<type_wrapper<vector<string> > >("string_array")
                .def(constructor<>())
        ]
    ];
    return 0;
}

} // namespace halmd
