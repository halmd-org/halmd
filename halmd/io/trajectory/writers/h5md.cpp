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
#include <halmd/io/trajectory/writers/h5md.hpp>
#include <halmd/io/utility/hdf5.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/demangle.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace io { namespace trajectory { namespace writers
{

/**
 * read sample from H5MD trajectory file
 */
template <int dimension, typename float_type>
h5md<dimension, float_type>::h5md(
    shared_ptr<sample_type> sample
  , shared_ptr<clock_type> clock
  , string const& file_name
)
  : sample(sample)
  , clock(clock)
  , file_(file_name, H5F_ACC_TRUNC)
  , time_(-1)
{
    LOG("write trajectory to file: " << file_.getFileName());

    // store file version in parameter group
    array<unsigned char, 2> version = {{ 1, 0 }};
    h5xx::write_attribute(h5xx::open_group(file_, "param"), "file_version", version);

    // open or create trajectory group
    H5::Group root = h5xx::open_group(file_, "trajectory");

    for (size_t i = 0; i < sample->r.size(); ++i) {
        H5::Group type;
        if (sample->r.size() > 1) {
            type = root.createGroup(string(1, 'A' + i));
        }
        else {
            type = root;
        }
        size_t size = sample->r[i]->size();
        // create and store datasets
        position_.push_back(h5xx::create_dataset<sample_vector_type>(type, "position", size));
        velocity_.push_back(h5xx::create_dataset<sample_vector_type>(type, "velocity", size));
    }

    // store data set for simulation time in reduced units
    time_ = h5xx::create_dataset<double>(root, "time");
}

/**
 * append samples to datasets
 */
template <int dimension, typename float_type>
void h5md<dimension, float_type>::append(uint64_t step)
{
    on_append_(step);

    if (sample->step != step) {
        throw logic_error("phase space sample was not updated");
    }

    // write particle positions
    for (size_t i = 0; i < position_.size(); ++i) {
        h5xx::write_dataset(position_[i], *sample->r[i]);
    }

    // write particle velocities
    for (size_t i = 0; i < velocity_.size(); ++i) {
        h5xx::write_dataset(velocity_[i], *sample->v[i]);
    }

    // write current simulation time
    h5xx::write_dataset(time_, clock->time());
}

/**
 * flush H5MD file
 */
template <int dimension, typename float_type>
void h5md<dimension, float_type>::flush()
{
    file_.flush(H5F_SCOPE_GLOBAL);
}

template <int dimension, typename float_type>
void h5md<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("h5md_" + lexical_cast<string>(dimension) + "_" + demangled_name<float_type>() + "_");
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("trajectory")
            [
                namespace_("writers")
                [
                    class_<h5md, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                            shared_ptr<sample_type>
                          , shared_ptr<clock_type>
                          , string const&
                        >())
                        .def("file", &h5md::file)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_trajectory_writers_h5md(lua_State* L)
{
    h5md<3, double>::luaopen(L);
    h5md<2, double>::luaopen(L);
    h5md<3, float>::luaopen(L);
    h5md<2, float>::luaopen(L);
    return 0;
}

}}} // namespace io::trajectory::writers

} // namespace halmd
