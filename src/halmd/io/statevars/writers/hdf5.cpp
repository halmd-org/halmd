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
#include <boost/filesystem/operations.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/io/statevars/writers/hdf5.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace std;
using namespace H5;

namespace halmd
{
namespace io { namespace statevars { namespace writers
{

/**
 * open HDF5 file for writing
 */
template <int dimension>
hdf5<dimension>::hdf5(modules::factory& factory, po::options const& vm)
  : _Base(factory, vm)
  , file_(
        (initial_path() / (vm["output"].as<string>() + ".msv")).string()
      , H5F_ACC_TRUNC // truncate existing file
    )
{
    // create parameter group
    Group param = open_group(file_, "/param");

    // store file version
    array<unsigned char, 2> version = {{ 1, 0 }};
    attribute(param, "file_version") = version;

    // store dimension
    // FIXME configuration and simulation parameters should be stored by a distinct function
    Group mdsim = param.createGroup("mdsim");
    attribute(mdsim, "dimension") = dimension;

    LOG("write macroscopic state variables to file: " << file_.getFileName());
}

/**
 * create dataset for scalar macroscopic state variable
 * and register writer functor
 */
template <int dimension> template<typename T>
void hdf5<dimension>::register_observable(
    string const& tag
  , T const* value_ptr
  , string const& desc
)
{
    // first part of tag is path, last part is dataset name
    list<string> path(split_path(tag));
    Group root = open_group(file_, path.begin(), --path.end());

    // create dataset for an unlimited number of scalar chunks
    DataSet dataset = create_dataset<T>(root, path.back());

    // store description as attribute
    attribute(dataset, "description") = desc;

    // add dataset writer to internal list
    writer_.push_back(make_dataset_writer(dataset, value_ptr));
}

template <int dimension> template<typename T>
void hdf5<dimension>::register_observable(
    string const& tag
  , vector<T> const* value_ptr
  , string const& desc
)
{
    // first part of tag is path, last part is dataset name
    list<string> path(split_path(tag));
    Group root = open_group(file_, path.begin(), --path.end());

    // create dataset for an unlimited number of vector chunks with given size
    DataSet dataset = create_dataset<vector<T> >(root, path.back(), value_ptr->size());

    // store description as attribute
    attribute(dataset, "description") = desc;

    // add dataset writer to internal list
    writer_.push_back(make_dataset_writer(dataset, value_ptr));
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

    SELECT_REGISTER_OBSERVABLE(double);
    SELECT_REGISTER_OBSERVABLE(vector_type);
    SELECT_REGISTER_OBSERVABLE(vector<double>);
    SELECT_REGISTER_OBSERVABLE(vector<vector_type>);
#undef SELECT_REGISTER_OBSERVABLE

    throw runtime_error(string("HDF5 writer: unknown type of dataset ") + tag);
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

}}} // namespace io::profile::writers

template class module<io::statevars::writers::hdf5<2> >;
template class module<io::statevars::writers::hdf5<3> >;

} // namespace halmd
