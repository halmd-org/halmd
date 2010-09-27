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

using namespace boost;
using namespace std;
using namespace H5;

namespace halmd
{
namespace io { namespace trajectory { namespace readers
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::depends()
{
    modules::depends<_Self, sample_type>::required();
}

template <int dimension, typename float_type>
void hdf5<dimension, float_type>::select(po::variables_map const& vm)
{
    if (!H5File::isHdf5(vm["trajectory-file"].as<string>())) {
        throw unsuitable_module("not an HDF5 file: " + vm["trajectory-file"].as<string>());
    }
}

/**
 * read sample from HDF5 trajectory file
 */
template <int dimension, typename float_type>
hdf5<dimension, float_type>::hdf5(modules::factory& factory, po::variables_map const& vm)
  : _Base(factory, vm)
  // dependency injection
  , sample(modules::fetch<sample_type>(factory, vm))
{
    LOG("read trajectory file: " << path_);

    H5File file(path_, H5F_ACC_RDONLY);
    Group root = open_group(file, "trajectory");

    for (size_t i = 0; i < sample->r.size(); ++i) {
        Group type;
        if (sample->r.size() > 1) {
            type = open_group(root, string(1, 'A' + i));
        }
        else {
            type = root;
        }

        DataSet r;
        try {
            // backwards compatibility with r:R:v:t format
            //   r = reduced single- or double-precision positions,
            //   R = extended single- or double-precision positions,
            //   v = single- or double-precision velocities
            //   t = single- or double-precision simulation time
            r = type.openDataSet("R");
            LOG_WARNING("detected obsolete trajectory file format");
        }
        catch (GroupIException const& e)
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
        DataSet v = type.openDataSet("v");

        H5::read(r, &*sample->r[i], offset_);
        H5::read(v, &*sample->v[i], offset_);
    }

    DataSet t = root.openDataSet("t");
    float_type time;
    size_t offset = H5::read(t, &time, offset_);
    LOG("read trajectory sample at offset " << offset << " with t = " << time);
}

// explicit instantiation
template class hdf5<3, double>;
template class hdf5<3, float>;
template class hdf5<2, double>;
template class hdf5<2, float>;

}}} // namespace io::trajectory::readers

template class module<io::trajectory::readers::hdf5<3, double> >;
template class module<io::trajectory::readers::hdf5<3, float> >;
template class module<io::trajectory::readers::hdf5<2, double> >;
template class module<io::trajectory::readers::hdf5<2, float> >;

} // namespace halmd
