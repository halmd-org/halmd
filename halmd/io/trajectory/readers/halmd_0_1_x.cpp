/*
 * Copyright © 2008-2011  Peter Colberg
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

#include <halmd/io/logger.hpp>
#include <halmd/io/trajectory/readers/halmd_0_1_x.hpp>
#include <halmd/io/utility/hdf5.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace io { namespace trajectory { namespace readers
{

/**
 * check whether file is in HALMD 0.1.x format
 *
 * FIXME implement file format detection based on required datasets
 */
template <int dimension, typename float_type>
optional<H5::H5File> halmd_0_1_x<dimension, float_type>::format(string const& file_name)
{
    if (!H5::H5File::isHdf5(file_name)) {
        return optional<H5::H5File>();
    }

    H5::H5File file(file_name, H5F_ACC_RDONLY);
    H5::Group param(file.openGroup("param"));

    if (h5xx::exists_attribute(param, "file_version")) {
        return optional<H5::H5File>();
    }

    LOG_DEBUG("detected HALMD 0.1.x trajectory file format");

    return file;
}

/**
 * read sample from HALMD 0.1.x trajectory file
 */
template <int dimension, typename float_type>
halmd_0_1_x<dimension, float_type>::halmd_0_1_x(
    shared_ptr<sample_type> sample
  , string const& file_name
  , ssize_t offset
)
  : sample(sample)
  , path_(file_name)
  , offset_(offset)
{
    LOG("read trajectory file: " << path_);

    H5::H5File file(path_, H5F_ACC_RDONLY);
    H5::Group root = h5xx::open_group(file, "trajectory");

    for (size_t i = 0; i < sample->r.size(); ++i) {
        H5::Group type;
        if (sample->r.size() > 1) {
            type = h5xx::open_group(root, string(1, 'A' + i));
        }
        else {
            type = root;
        }

        H5::DataSet r;
        try {
            // backwards compatibility with r:R:v:t format
            //   r = reduced single- or double-precision positions,
            //   R = extended single- or double-precision positions,
            //   v = single- or double-precision velocities
            //   t = single- or double-precision simulation time
            r = type.openDataSet("R");
            LOG_WARNING("detected obsolete trajectory file format");
        }
        catch (H5::GroupIException const& e)
        {
            // new-style r:v:t format
            //   r = extended double-precision positions,
            //   v = single- or double-precision velocities
            //   t = double-precision simulation time
            r = type.openDataSet("r");
        }
        // backwards compatibility with r:R:v:t format
        if (r.getDataType().getId() == h5xx::ctype<float>::hid()) {
            // use reduced positions if extended positions are single-precision
            r = type.openDataSet("r");
            LOG_WARNING("falling back to reduced particle position sample");
        }
        H5::DataSet v = type.openDataSet("v");

        h5xx::read_dataset(r, sample->r[i].get(), offset_);
        h5xx::read_dataset(v, sample->v[i].get(), offset_);
    }

    H5::DataSet t = root.openDataSet("t");
    float_type time;
    offset = h5xx::read_dataset(t, &time, offset_);
    LOG("read trajectory sample at offset " << offset << " with t = " << time);
}

template <int dimension, typename float_type>
void halmd_0_1_x<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("halmd_0_1_x_" + lexical_cast<string>(dimension) + "_" + demangled_name<float_type>() + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("io")
            [
                namespace_("trajectory")
                [
                    namespace_("readers")
                    [
                        class_<halmd_0_1_x, shared_ptr<_Base>, _Base>(class_name.c_str())
                            .def(constructor<
                                shared_ptr<sample_type>
                                , string const&
                                , ssize_t
                            >())
                            .scope
                            [
                                def("format", &halmd_0_1_x::format)
                            ]
                    ]
                ]
            ]
        ]
    ];
}

HALMD_INIT( register_luaopen )
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &halmd_0_1_x<3, double>::luaopen
    ]
    [
        &halmd_0_1_x<2, double>::luaopen
    ]
    [
        &halmd_0_1_x<3, float>::luaopen
    ]
    [
        &halmd_0_1_x<2, float>::luaopen
    ];
}

}}} // namespace io::trajectory::readers

} // namespace halmd
