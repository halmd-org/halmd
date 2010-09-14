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

#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/operations.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/io/profile/writers/hdf5.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace std;

namespace halmd
{
namespace io { namespace profile { namespace writers
{

/**
 * open HDF5 file for writing
 */
hdf5::hdf5(modules::factory& factory, po::options const& vm)
  : _Base(factory, vm)
  , file_(
        (initial_path() / (vm["output"].as<string>() + ".prf")).string()
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

    LOG("write profile data to file: " << file_.getFileName());
}

/**
 * create dataset for runtime accumulator
 */
void hdf5::register_accumulator(
    vector<string> const& tag
  , accumulator_type const& acc
  , string const& desc
)
{
    H5::Group group = file_.openGroup("/");
    vector<string>::const_iterator it = tag.begin();
    vector<string>::const_iterator const end = tag.end() - 1;
    while (it != end) { // create group for each tag token
        try {
            H5XX_NO_AUTO_PRINT(H5::GroupIException);
            group = group.openGroup(*it);
        }
        catch (H5::GroupIException const&) {
            group = group.createGroup(*it);
        }
        ++it;
    }
    array<hsize_t, 1> dim = {{ 3 }};
    H5::DataSet dataset(
        group.createDataSet(
            trim_right_copy_if(tag.back(), is_any_of("_")) // omit trailing "_"
          , H5::PredType::NATIVE_DOUBLE
          , H5::DataSpace(dim.size(), dim.data())
        )
    );
    H5::Attribute attr( // store description as attribute
        dataset.createAttribute(
            "timer"
          , H5::StrType(H5::PredType::C_S1, desc.size()) // no NULL terminator
          , H5S_SCALAR
        )
    );
    attr.write(attr.getDataType(), desc.c_str());

    // We bind the functions to write the datasets using a
    // *reference* to the accumulator and a *copy* of the HDF5
    // dataset instance which goes out of scope

    writer_.push_back(
        bind(
            &hdf5::write_accumulator
          , dataset
          , cref(acc)
        )
    );
}

/**
 * write dataset for runtime accumulator
 */
void hdf5::write_accumulator(
    H5::DataSet const& dataset
  , accumulator_type const& acc
)
{
    array<hsize_t, 1> dim = {{ 3 }};
    array<double, 3> data = {{
        mean(acc)
      , error_of_mean(acc)
      , count(acc)
    }};
    dataset.write(
        data.data()
      , H5::PredType::NATIVE_DOUBLE
      , H5::DataSpace(dim.size(), dim.data())
      , dataset.getSpace()
    );
}

/**
 * write all datasets and flush file to disk
 */
void hdf5::write()
{
    for_each(
        writer_.begin()
      , writer_.end()
      , bind(&writer_functor::operator(), _1)
    );
    file_.flush(H5F_SCOPE_GLOBAL);
}

}}} // namespace io::profile::writers

template class module<io::profile::writers::hdf5>;

} // namespace halmd
