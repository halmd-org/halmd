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
    H5::Group param = H5::open_group(file_, "/param");

    // store file version
    array<unsigned char, 2> version = {{ 1, 0 }};
    H5::attribute(param, "file_version") = version;

    // store dimension
    // FIXME configuration and simulation parameters should be stored by a distinct function
    H5::Group mdsim = param.createGroup("mdsim");
    H5::attribute(mdsim, "dimension") = dimension;

    LOG("write macroscopic state variables to file: " << file_.getFileName());
}

/**
 * create dataset for scalar macroscopic state variable
 * and register writer functor
 */
template <int dimension>
void hdf5<dimension>::register_scalar_observable(
    string const& tag
  , double const* value_ptr
  , string const& desc
)
{
    H5::Group root = H5::open_group(file_, "/");

    // create dataset for an unlimited number of scalar chunks
    H5::DataSet dataset = H5::create_dataset<double>(root, tag);

    // store description as attribute
    H5::attribute(dataset, "description") = desc;

    // add dataset writer to internal list
    writer_.push_back(H5::make_dataset_writer(dataset, value_ptr));
}

/**
 * create dataset for vectorial macroscopic state variable
 * and register writer functor
 */
template <int dimension>
void hdf5<dimension>::register_vector_observable(
    string const& tag
  , vector_type const* value_ptr
  , string const& desc
)
{
    // FIXME overloading with fixed_vector is missing in H5
    // shouldn't fixed_vector be interpreted as boost::array?
    typedef boost::array<typename vector_type::value_type, vector_type::static_size> array_type;
    H5::Group root = H5::open_group(file_, "/");

    // create dataset for an unlimited number of vectorial chunks
    H5::DataSet dataset = H5::create_dataset<array_type>(root, tag);

    // store description as attribute
    H5::attribute(dataset, "description") = desc;

    // add dataset writer to internal list
    writer_.push_back(H5::make_dataset_writer(dataset, static_cast<array_type const*>(value_ptr)));
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
