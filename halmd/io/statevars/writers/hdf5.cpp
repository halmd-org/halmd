/*
 * Copyright © 2010  Peter Colberg and Felix Höfling
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

#include <boost/algorithm/string/trim.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/io/statevars/writers/hdf5.hpp>
#include <halmd/io/utility/hdf5.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/demangle.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace std;

namespace halmd
{
namespace io { namespace statevars { namespace writers
{

/**
 * open HDF5 file for writing
 */
template <int dimension>
hdf5<dimension>::hdf5(string const& file_name)
  : file_(file_name, H5F_ACC_TRUNC) // truncate existing file
{
    // create parameter group
    H5::Group param = h5xx::open_group(file_, "/param");

    // store file version
    array<unsigned char, 2> version = {{ 1, 0 }};
    h5xx::write_attribute(param, "file_version", version);

    LOG("write macroscopic state variables to file: " << file_.getFileName());
}

/**
 * create dataset for macroscopic state variable
 * and register writer functor
 */
template <int dimension> template<typename T>
void hdf5<dimension>::register_observable(
    string const& tag
  , T const* value_ptr
  , string const& desc
)
{
    // create dataset for an unlimited number of chunks
    // first part of tag is path, last part is dataset name
    H5::DataSet dataset = h5xx::create_dataset<T>(file_, tag);

    // store description as attribute
    h5xx::write_attribute(dataset, "description", desc);

    // add dataset writer to internal list
    writer_.push_back(h5xx::make_dataset_writer(dataset, value_ptr));
}

// overload for vector of data
template <int dimension> template<typename T>
void hdf5<dimension>::register_observable(
    string const& tag
  , vector<T> const* value_ptr
  , string const& desc
)
{
    // create dataset for an unlimited number of vector chunks with given size
    // first part of tag is path, last part is dataset name
    H5::DataSet dataset = h5xx::create_dataset<vector<T> >(file_, tag, value_ptr->size());

    // store description as attribute
    h5xx::write_attribute(dataset, "description", desc);

    // add dataset writer to internal list
    writer_.push_back(h5xx::make_dataset_writer(dataset, value_ptr));
}

/**
 * register function from base class interface,
 * calls appropriate template of register_observable
 */
template <int dimension>
void hdf5<dimension>::register_observable(
    string const& tag
  , void const* value_ptr
  , type_info const& value_type
  , string const& desc
)
{
    // the reinterpret_cast is safe since we know the type of value_ptr
#define SELECT_REGISTER_OBSERVABLE(TYPE)               \
    if (value_type == typeid(TYPE)) {                  \
        register_observable<>(                         \
            tag                                        \
          , reinterpret_cast<TYPE const*>(value_ptr)   \
          , desc                                       \
        );                                             \
        return;                                        \
    }

    typedef boost::array<double, 3> array_double_3_type; // avoid comma in macro parameter

    SELECT_REGISTER_OBSERVABLE(double);
    SELECT_REGISTER_OBSERVABLE(unsigned int);
    SELECT_REGISTER_OBSERVABLE(vector_type);
    SELECT_REGISTER_OBSERVABLE(array_double_3_type);
    SELECT_REGISTER_OBSERVABLE(vector<double>);
    SELECT_REGISTER_OBSERVABLE(vector<unsigned int>);
    SELECT_REGISTER_OBSERVABLE(vector<vector_type>);
    SELECT_REGISTER_OBSERVABLE(vector<array_double_3_type>);
#undef SELECT_REGISTER_OBSERVABLE

    throw runtime_error(
        string("HDF5 writer: unknown type of dataset ") + tag
      + " (" + demangled_name(value_type) + ")"
    );
}

/**
 * create single-valued dataset and write data
 */
template <int dimension> template<typename T>
void hdf5<dimension>::write_dataset(
    string const& tag
  , T const& value
  , string const& desc
)
{
    // create dataset for a single chunk
    // first part of tag is path, last part is dataset name
    H5::DataSet dataset = h5xx::create_dataset<T>(file_, tag, 1);

    // store description as attribute
    h5xx::write_attribute(dataset, "description", desc);

    // write dataset at index 0
    h5xx::write_dataset(dataset, value, 0);
}

// overload for vector of data
template <int dimension> template<typename T>
void hdf5<dimension>::write_dataset(
    string const& tag
  , vector<T> const& value
  , string const& desc
)
{
    // create dataset for a single vector chunk with given size
    // first part of tag is path, last part is dataset name
    H5::DataSet dataset = h5xx::create_dataset<vector<T> >(file_, tag, value.size(), 1);

    // store description as attribute
    h5xx::write_attribute(dataset, "description", desc);

    // write dataset at index 0
    h5xx::write_dataset(dataset, value, 0);
}

/**
 * write_dataset function from base class interface,
 * calls appropriate template of write_dataset
 */
template <int dimension>
void hdf5<dimension>::write_dataset(
    string const& tag
  , void const* value_ptr
  , type_info const& value_type
  , string const& desc
)
{
    // the reinterpret_cast is safe since we know the type of value_ptr
#define SELECT_WRITE_DATASET(TYPE)                     \
    if (value_type == typeid(TYPE)) {                  \
        write_dataset<>(                               \
            tag                                        \
          , *reinterpret_cast<TYPE const*>(value_ptr)  \
          , desc                                       \
        );                                             \
        return;                                        \
    }

    typedef boost::array<double, 3> array_double_3_type; // avoid comma in macro parameter

    SELECT_WRITE_DATASET(double);
    SELECT_WRITE_DATASET(unsigned int);
    SELECT_WRITE_DATASET(vector_type);
    SELECT_WRITE_DATASET(array_double_3_type);
    SELECT_WRITE_DATASET(vector<double>);
    SELECT_WRITE_DATASET(vector<unsigned int>);
    SELECT_WRITE_DATASET(vector<vector_type>);
    SELECT_WRITE_DATASET(vector<array_double_3_type>);
#undef SELECT_WRITE_DATASET

    throw runtime_error(
        string("HDF5 writer: unknown type of dataset ") + tag
      + " (" + demangled_name(value_type) + ")"
    );
}

/**
 * write all datasets and flush file to disk
 */
template <int dimension>
void hdf5<dimension>::write()
{
    for_each(
        writer_.begin()
      , writer_.end()
      , bind(&writer_functor::operator(), _1)
    );
    file_.flush(H5F_SCOPE_GLOBAL);
}

template <int dimension>
void hdf5<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("hdf5_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("statevars")
            [
                namespace_("writers")
                [
                    class_<hdf5, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<string const&>())
                        .def("file", &hdf5::file)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_statevars_writers_hdf5(lua_State* L)
{
    hdf5<3>::luaopen(L);
    hdf5<2>::luaopen(L);
    return 0;
}

}}} // namespace io::statevars::writers

} // namespace halmd