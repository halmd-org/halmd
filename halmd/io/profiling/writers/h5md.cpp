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

#include <halmd/io/logger.hpp>
#include <halmd/io/profiling/writers/h5md.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace io { namespace profiling { namespace writers
{

/**
 * open H5MD file for writing
 */
h5md::h5md(string const& file_name)
  : file_(file_name, H5F_ACC_TRUNC) // truncate existing file
{
    // create parameter group
    H5::Group param = h5xx::open_group(file_, "/param");

    // store file version
    array<unsigned char, 2> version = {{ 1, 0 }};
    h5xx::write_attribute(param, "file_version", version);

    LOG("write profiling data to file: " << file_.getFileName());
}

/**
 * create dataset for runtime accumulator
 */
void h5md::register_accumulator(
    string const& group
  , accumulator_type const& acc
  , string const& tag
  , string const& desc
)
{
    // open group defined by tag,
    // last entry of tag will be the name of the dataset
    H5::DataSet dataset = h5xx::create_dataset<array<double, 3> >(
        file_
      , "profiling/" + group + "/" + tag
      , 1                                                   // only 1 entry
    );
    // store description as attribute
    h5xx::write_attribute(dataset, "timer", desc);

    // We bind the functions to write the datasets using a
    // *reference* to the accumulator and a *copy* of the HDF5
    // dataset instance which goes out of scope

    writer_.push_back(
        bind(
            &h5md::write_accumulator
          , dataset
          , cref(acc)
        )
    );
}

/**
 * write dataset for runtime accumulator
 */
void h5md::write_accumulator(
    H5::DataSet const& dataset
  , accumulator_type const& acc
)
{
    array<double, 3> data = {{
        mean(acc)
      , error_of_mean(acc)
      , count(acc)
    }};
    h5xx::write_dataset(dataset, data, 0);
}

/**
 * write all datasets and flush file to disk
 */
void h5md::write()
{
    for_each(
        writer_.begin()
      , writer_.end()
      , bind(&writer_functor::operator(), _1)
    );
    file_.flush(H5F_SCOPE_GLOBAL);
}

void h5md::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("io")
            [
                namespace_("profiling")
                [
                    namespace_("writers")
                    [
                        class_<h5md, shared_ptr<_Base>, _Base>("h5md")
                            .def(constructor<string const&>())
                            .def("file", &h5md::file)
                    ]
                ]
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &h5md::luaopen
    ];
}

} // namespace

}}} // namespace io::profiling::writers

} // namespace halmd
