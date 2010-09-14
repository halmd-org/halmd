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
      , H5xx::file::trunc // truncate existing file
    )
{
    // create parameter group
    H5::Group param = file_.openGroup("/").createGroup("param");

    // store file version
    array<hsize_t, 1> dim = {{ 2 }};
    array<unsigned char, 2> version = {{ 1, 0 }};
    H5::Attribute attr(
        param.createAttribute(
            "file_version"
          , H5::PredType::NATIVE_UCHAR
          , H5::DataSpace(dim.size(), dim.data())
        )
    );
    attr.write(attr.getDataType(), version.data());

    // store dimension
    // FIXME configuration and simulation parameters should be stored by a distinct function
    array<unsigned, 1> dimension_ = {{ dimension }};
    H5::Group mdsim = param.createGroup("mdsim");
    attr = mdsim.createAttribute("dimension", H5::PredType::NATIVE_UINT, H5S_SCALAR);
    attr.write(attr.getDataType(), dimension_.data());

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
    H5::Group root = file_.openGroup("/");

    // write an unlimited number of scalar chunks
    array<hsize_t, 1> dim = {{ 0 }};
    array<hsize_t, 1> max_dim = {{ H5S_UNLIMITED }};
    array<hsize_t, 1> chunk_dim = {{ 1 }};

    H5::DSetCreatPropList cparms;
    cparms.setChunk(chunk_dim.size(), chunk_dim.data());

    H5::DataSet dataset(
        root.createDataSet(
            tag
          , H5::PredType::NATIVE_DOUBLE
          , H5::DataSpace(dim.size(), dim.data(), max_dim.data())
          , cparms
        )
    );

    // store description as attribute
    H5::Attribute attr(
        dataset.createAttribute(
            "description"
          , H5::StrType(H5::PredType::C_S1, desc.size()) // no NULL terminator
          , H5S_SCALAR
        )
    );
    attr.write(attr.getDataType(), desc.c_str());

    // We bind the functions to write the datasets, using a
    // *pointer* to the value and a *copy* of the HDF5
    // dataset instance which goes out of scope
    writer_.push_back(bind(&hdf5<dimension>::write_scalar_observable, dataset, value_ptr));
}

/**
 * write dataset for scalar macroscopic state variable
 */
template <int dimension>
void hdf5<dimension>::write_scalar_observable(
    H5::DataSet const& dataset
  , double const* value_ptr
)
{
    array<hsize_t, 1> dim;

    // extend data space by 1
    H5::DataSpace dataspace(dataset.getSpace());
    dataspace.getSimpleExtentDims(dim.data());
    hsize_t count[1]  = { 1 };
    hsize_t start[1]  = { dim[0] };
    hsize_t stride[1] = { 1 };
    hsize_t block[1]  = { 1 };
    dim[0]++;
    dataspace.setExtentSimple(dim.size(), dim.data());
    dataspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    dataset.extend(dim.data());
    dataset.write(value_ptr, H5::PredType::NATIVE_DOUBLE, H5S_SCALAR, dataspace);
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
    H5::Group root = file_.openGroup("/");

    // write an unlimited number of vectorial chunks
    array<hsize_t, 2> dim = {{ 0, vector_type::static_size }};
    array<hsize_t, 2> max_dim = {{ H5S_UNLIMITED, vector_type::static_size }};
    array<hsize_t, 2> chunk_dim = {{ 1, vector_type::static_size }};

    H5::DSetCreatPropList cparms;
    cparms.setChunk(chunk_dim.size(), chunk_dim.data());

    H5::DataSet dataset(
        root.createDataSet(
            tag
          , H5xx::ctype<typename vector_type::value_type>()
          , H5::DataSpace(dim.size(), dim.data(), max_dim.data())
          , cparms
        )
    );

    // store description as attribute
    H5::Attribute attr(
        dataset.createAttribute(
            "description"
          , H5::StrType(H5::PredType::C_S1, desc.size()) // no NULL terminator
          , H5S_SCALAR
        )
    );
    attr.write(attr.getDataType(), desc.c_str());

    // We bind the functions to write the datasets, using a
    // *pointer* to the value and a *copy* of the HDF5
    // dataset instance which goes out of scope
    writer_.push_back(bind(&hdf5<dimension>::write_vector_observable, dataset, value_ptr));
}

/**
 * write dataset for scalar macroscopic state variable
 */
template <int dimension>
void hdf5<dimension>::write_vector_observable(
    H5::DataSet const& dataset
  , vector_type const* value_ptr
)
{
    array<hsize_t, 2> dim;

    // extend data space by 1
    H5::DataSpace dataspace(dataset.getSpace());
    dataspace.getSimpleExtentDims(dim.data());
    hsize_t count[2]  = { 1, 1 };
    hsize_t start[2]  = { dim[0], 0 };
    hsize_t stride[2] = { 1, 1 };
    hsize_t block[2]  = { 1, vector_type::static_size };
    dim[0]++;
    dataspace.setExtentSimple(dim.size(), dim.data());
    dataspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    array<hsize_t, 1> dim_sample = {{ vector_type::static_size }};

    dataset.extend(dim.data());
    dataset.write(
        value_ptr->data()
      , H5xx::ctype<typename vector_type::value_type>()
      , H5::DataSpace(dim_sample.size(), dim_sample.data())
      , dataspace
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

}}} // namespace io::profile::writers

template class module<io::statevars::writers::hdf5<2> >;
template class module<io::statevars::writers::hdf5<3> >;

} // namespace halmd
