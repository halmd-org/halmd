/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <boost/foreach.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/io/trajectory/writers/hdf5.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace boost::filesystem;
using namespace std;

namespace halmd
{
namespace io { namespace trajectory { namespace writers
{

/**
 * read sample from HDF5 trajectory file
 */
template <int dimension, typename float_type>
hdf5<dimension, float_type>::hdf5(
    shared_ptr<sample_type> sample
  , string const& file_name
)
  : sample(sample)
  , path_(initial_path() / file_name)
  , file_(path_.file_string(), H5F_ACC_TRUNC)
{
    LOG("write trajectory to file: " << path_.file_string());

    // store file version in parameter group
    array<unsigned char, 2> version = {{ 1, 0 }};
    H5::attribute(open_group(file_, "param"), "file_version") = version;

    // open or create trajectory group
    H5::Group root = open_group(file_, "trajectory");

    for (size_t i = 0; i < sample->r.size(); ++i) {
        H5::Group type;
        if (sample->r.size() > 1) {
            type = root.createGroup(string(1, 'A' + i));
        }
        else {
            type = root;
        }
        size_t size = sample->r[i]->size();
        H5::DataSet r = H5::create_dataset<sample_vector_type>(type, "position", size);
        H5::DataSet v = H5::create_dataset<sample_vector_type>(type, "velocity", size);

        // particle positions
        writers_.push_back(make_dataset_writer(r, &*sample->r[i]));
        // particle velocities
        writers_.push_back(make_dataset_writer(v, &*sample->v[i]));
    }

    // simulation time in reduced units
    H5::DataSet t = H5::create_dataset<double>(root, "time");
    writers_.push_back(make_dataset_writer(t, &sample->time));
}

/**
 * append samples to datasets
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::append()
{
    sample->acquire();
    BOOST_FOREACH (boost::function<void ()> const& writer, writers_) {
        writer();
    }
}

/**
 * flush HDF5 file
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::flush()
{
    file_.flush(H5F_SCOPE_GLOBAL);
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
                    namespace_("writers")
                    [
                        class_<T, shared_ptr<_Base>, _Base>(class_name)
                            .def(constructor<
                                shared_ptr<sample_type>
                              , string const&
                            >())
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

}}} // namespace io::trajectory::writers

} // namespace halmd
