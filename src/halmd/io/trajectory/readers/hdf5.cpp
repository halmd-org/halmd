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

#include <halmd/io/logger.hpp>
#include <halmd/io/trajectory/readers/hdf5.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace io { namespace trajectory { namespace readers
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void hdf5<dimension, float_type>::resolve(po::options const& vm)
{
    module<sample_type>::required(vm);
    module<particle_type>::required(vm);
}

/**
 * read sample from HDF5 trajectory file
 */
template <int dimension, typename float_type>
hdf5<dimension, float_type>::hdf5(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , sample(module<sample_type>::fetch(vm))
  , particle(module<particle_type>::fetch(vm))
  // parse options
  , path_(vm["trajectory"].as<string>())
  , offset_(vm["trajectory-sample"].as<ssize_t>())
{
    LOG("read trajectory file: " << path_);

    H5::H5File file(path_, H5F_ACC_RDONLY);
    H5::Group root(file.openGroup("trajectory"));

    for (size_t i = 0; i < particle->ntype; ++i) {
        H5::Group type;
        if (particle->ntype > 1) {
            type = root.openGroup(string(1, 'A' + i));
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
        if (r.getDataType() == H5xx::ctype<float>()) {
            // use reduced positions if extended positions are single-precision
            r = type.openDataSet("r");
            LOG_WARNING("falling back to reduced particle position sample");
        }
        H5::DataSet v = type.openDataSet("v");

        sample->r[i].reset(new position_sample_vector(particle->ntypes[i]));
        sample->v[i].reset(new velocity_sample_vector(particle->ntypes[i]));

        read(r, sample->r[i]);
        read(v, sample->v[i]);
    }

    H5::DataSet t = root.openDataSet("t");
    float_type time;
    size_t offset = read(t, time);
    LOG("read trajectory sample at offset " << offset << " with t = " << time);
}

/**
 * read vector sample dataset
 */
template <int dimension, typename float_type>
size_t hdf5<dimension, float_type>::read(H5::DataSet dset, position_sample_ptr sample)
{
    H5::DataSpace ds(dset.getSpace());

    if (!ds.isSimple()) {
        throw runtime_error("HDF5 vector dataspace is not a simple dataspace");
    }
    if (ds.getSimpleExtentNdims() != 3) {
        throw runtime_error("HDF5 vector dataspace has invalid dimensionality");
    }

    hsize_t dim[3];
    ds.getSimpleExtentDims(dim);

    ssize_t const len = dim[0];
    if ((offset_ >= len) || ((-offset_) > len)) {
        throw runtime_error("trajectory input sample number out of bounds");
    }
    size_t offset = (offset_ < 0) ? (offset_ + len) : offset_;

    size_t const size = dim[1];
    if (size != sample->size()) {
        throw runtime_error("trajectory input file has invalid number of particles");
    }
    if (dim[2] != dimension) {
        throw runtime_error("trajectory input file has invalid coordinate dimension");
    }

    hsize_t dim_sample[2] = { size, dimension };
    H5::DataSpace ds_sample(2, dim_sample);

    hsize_t count[3]  = { 1, size, 1 };
    hsize_t start[3]  = { offset, 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, 1, dimension };
    ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    try {
        H5XX_NO_AUTO_PRINT(H5::Exception);
        dset.read(sample->data(), H5xx::ctype<float_type>(), ds_sample, ds);
    }
    catch (H5::Exception const&) {
        throw runtime_error("failed to read vector sample from HDF5 trajectory input file");
    }

    return offset;
}

/**
 * read scalar sample dataset
 */
template <int dimension, typename float_type>
size_t hdf5<dimension, float_type>::read(H5::DataSet dset, float_type& sample)
{
    H5::DataSpace ds(dset.getSpace());

    if (!ds.isSimple()) {
        throw runtime_error("HDF5 scalar dataspace is not a simple dataspace");
    }
    if (ds.getSimpleExtentNdims() != 1) {
        throw runtime_error("HDF5 scalar dataspace has invalid dimensionality");
    }

    hsize_t dim[1];
    ds.getSimpleExtentDims(dim);

    ssize_t const len = dim[0];
    if ((offset_ >= len) || ((-offset_) > len)) {
        throw runtime_error("trajectory input sample number out of bounds");
    }
    size_t offset = (offset_ < 0) ? (offset_ + len) : offset_;

    hsize_t count[1]  = { 1 };
    hsize_t start[1]  = { offset };
    hsize_t stride[1] = { 1 };
    hsize_t block[1]  = { 1 };
    ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    try {
        H5XX_NO_AUTO_PRINT(H5::Exception);
        dset.read(&sample, H5xx::ctype<float_type>(), H5S_SCALAR, ds);
    }
    catch (H5::Exception const&) {
        throw runtime_error("failed to read scalar sample from HDF5 trajectory input file");
    }

    return offset;
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
