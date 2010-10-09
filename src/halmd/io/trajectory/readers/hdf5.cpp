/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <H5xx.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/io/trajectory/readers/hdf5.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace io { namespace trajectory { namespace readers
{

template <int dimension, typename float_type>
bool hdf5<dimension, float_type>::check(std::string const& file_name)
{
    return H5::H5File::isHdf5(file_name);
}

/**
 * read sample from HDF5 trajectory file
 */
template <int dimension, typename float_type>
hdf5<dimension, float_type>::hdf5(
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
    H5::Group root = open_group(file, "trajectory");

    for (size_t i = 0; i < sample->r.size(); ++i) {
        H5::Group type;
        if (sample->r.size() > 1) {
            type = open_group(root, string(1, 'A' + i));
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
        if (r.getDataType() == H5::ctype<float>()) {
            // use reduced positions if extended positions are single-precision
            r = type.openDataSet("r");
            LOG_WARNING("falling back to reduced particle position sample");
        }
        H5::DataSet v = type.openDataSet("v");

        H5::read(r, &*sample->r[i], offset_);
        H5::read(v, &*sample->v[i], offset_);
    }

    H5::DataSet t = root.openDataSet("t");
    float_type time;
    offset = H5::read(t, &time, offset_);
    LOG("read trajectory sample at offset " << offset << " with t = " << time);
}

template <typename T>
static void register_lua(lua_State* L, char const* class_name)
{
    typedef typename T::_Base _Base;
    typedef typename T::sample_type sample_type;

    using namespace luabind;
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
                        class_<T, shared_ptr<_Base>, _Base>(class_name)
                            .def(constructor<
                                shared_ptr<sample_type>
                                , string const&
                                , ssize_t
                            >())
                            .scope
                            [
                                def("check", &T::check)
                            ]
                    ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        bind(&register_lua<hdf5<3, double> >, _1, "hdf5_3_double_")
    ]
    [
        bind(&register_lua<hdf5<2, double> >, _1, "hdf5_2_double_")
    ]
    [
        bind(&register_lua<hdf5<3, float> >, _1, "hdf5_3_float_")
    ]
    [
        bind(&register_lua<hdf5<2, float> >, _1, "hdf5_2_float_")
    ];
}

}}} // namespace io::trajectory::readers

} // namespace halmd
