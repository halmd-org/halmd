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
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace std;
using namespace H5;

namespace halmd
{
namespace io { namespace profile { namespace writers
{

/**
 * open HDF5 file for writing
 */
hdf5::hdf5(string const& file_name)
  : file_(
        (initial_path() / file_name).string()
      , H5F_ACC_TRUNC // truncate existing file
    )
{
    // create parameter group
    Group param = open_group(file_, "/param");

    // store file version
    array<unsigned char, 2> version = {{ 1, 0 }};
    attribute(param, "file_version") = version;

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
    // open group defined by tag,
    // last entry of tag will be the name of the attribute
    Group group = open_group(file_, tag.begin(), tag.end() - 1);

    DataSet dataset = create_dataset<array<double, 3> >(
        group
      , trim_right_copy_if(tag.back(), is_any_of("_"))      // omit trailing "_"
      , 1                                                   // only 1 entry
    );
    // store description as attribute
    attribute(dataset, "timer") = desc;

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
    DataSet const& dataset
  , accumulator_type const& acc
)
{
    array<double, 3> data = {{
        mean(acc)
      , error_of_mean(acc)
      , count(acc)
    }};
    H5::write(dataset, data, 0);
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

static __attribute__((constructor)) void register_lua()
{
    typedef hdf5::_Base _Base;

    using namespace luabind;
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        namespace_("halmd_wrapper")
        [
            namespace_("io")
            [
                namespace_("profile")
                [
                    namespace_("writers")
                    [
                        class_<hdf5, shared_ptr<hdf5>, _Base>("hdf5")
                            .def(constructor<string const&>())
                    ]
                ]
            ]
        ]
    ];
}

}}} // namespace io::profile::writers

} // namespace halmd
