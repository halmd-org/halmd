/* MD simulation trajectory writer
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef LJGPU_SAMPLE_TRAJECTORY_HPP
#define LJGPU_SAMPLE_TRAJECTORY_HPP

#include <H5Cpp.h>
#include <algorithm>
#include <assert.h>
#include <boost/array.hpp>
#include <boost/type_traits.hpp>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/mdsim/sample.hpp>
#include <ljgpu/util/H5param.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <ljgpu/util/exception.hpp>
#include <ljgpu/util/log.hpp>
#include <string>

namespace ljgpu {

/**
 * trajectory file writer
 */
class trajectory
{
public:
    /** io flags */
    enum openmode {
	in = 0x1,
	out = 0x2,
    };

public:
    /** create HDF5 trajectory output file */
    void open(std::string const& filename, openmode mode = in);
    /** close HDF5 trajectory output file */
    void close();
    /** flush HDF5 output file to disk */
    void flush();

    /** read phase space sample */
    template <typename sample_type>
    void read(sample_type& sample, int64_t index);
    /** write phase space sample */
    template <typename sample_type>
    void write(sample_type const& sample, double time);

    /** returns HDF5 parameter group */
    H5param attrs();

private:
    /** HDF5 trajectory output file */
    H5::H5File m_file;
};

/**
 * read phase space sample
 */
template <typename sample_type>
void trajectory::read(sample_type& sample, int64_t index)
{
    enum { dimension = sample_type::position_vector::static_size };
    typedef typename sample_type::position_vector::value_type position_value_type;
    typedef typename sample_type::velocity_vector::value_type velocity_value_type;

    H5::DataType const tid_r(H5xx::ctype<position_value_type>::type);
    H5::DataType const tid_v(H5xx::ctype<velocity_value_type>::type);

    H5::DataSet dset_r, dset_v;
    try {
	H5::Group root(m_file.openGroup("trajectory"));
	try {
	    // backwards compatibility with r:R:v:t format
	    // 	 r = reduced single- or double-precision positions,
	    //   R = extended single- or double-precision positions,
	    //   v = single- or double-precision velocities
	    //   t = single- or double-precision simulation time
	    dset_r = root.openDataSet("R");
	    LOG_WARNING("detected obsolete trajectory file format");
	}
	catch (H5::Exception const&)
	{
	    // new-style r:v:t format
	    //   r = extended double-precision positions,
	    //   v = single- or double-precision velocities
	    //   t = double-precision simulation time
	    dset_r = root.openDataSet("r");
	}
	// backwards compatibility with r:R:v:t format
	if (dset_r.getDataType() == H5::PredType::NATIVE_FLOAT) {
	    // use reduced positions if extended positions are single-precision
	    dset_r = root.openDataSet("r");
	    LOG_WARNING("falling back to reduced particle position sample");
	}
	dset_v = root.openDataSet("v");
    }
    catch (H5::Exception const&) {
	throw exception("failed to open HDF5 trajectory datasets");
    }

    // validate dataspace extents
    H5::DataSpace ds_r(dset_r.getSpace());
    H5::DataSpace ds_v(dset_v.getSpace());
    if (!ds_r.isSimple()) {
	throw exception("trajectory dataspace is not a simple dataspace");
    }
    if (!ds_v.isSimple()) {
	throw exception("velocity dataspace is not a simple dataspace");
    }
    if (ds_r.getSimpleExtentNdims() != 3) {
	throw exception("trajectory dataspace has invalid dimensionality");
    }
    if (ds_v.getSimpleExtentNdims() != 3) {
	throw exception("velocity dataspace has invalid dimensionality");
    }

    // validate dataspace dimensions
    hsize_t dim_r[3];
    ds_r.getSimpleExtentDims(dim_r);
    hsize_t dim_v[3];
    ds_v.getSimpleExtentDims(dim_v);
    if (!std::equal(dim_r, dim_r + 3, dim_v)) {
	throw exception("trajectory and velocity dataspace dimensions differ");
    }

    // number of phase space samples
    int64_t len = dim_r[0];
    // number of particles
    unsigned int npart = dim_r[1];

    // validate dataspace dimensions
    if (len < 1) {
	throw exception("trajectory input file has invalid number of samples");
    }
    if (npart < 1) {
	throw exception("trajectory input file has invalid number of particles");
    }
    if (dimension != dim_r[2]) {
	throw exception("trajectory input file has invalid coordinate dimension");
    }

    // check if sample number is within bounds
    if ((index >= len) || ((-index) > len)) {
	throw exception("trajectory input sample number out of bounds");
    }
    index = (index < 0) ? (index + len) : index;

    LOG("resuming from trajectory sample at offset: " << index);

    try {
	sample.r.resize(npart);
	sample.v.resize(npart);
    }
    catch (std::bad_alloc const&) {
	throw exception("failed to allocate memory for trajectory input sample");
    }

    // read sample from dataset
    hsize_t dim_mem[2] = { npart, dimension };
    H5::DataSpace ds_mem(2, dim_mem);

    hsize_t count[3]  = { 1, npart, 1 };
    hsize_t start[3]  = { index, 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, 1, dimension };

    try {
	// read periodically reduced particle positions
	ds_r.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dset_r.read(sample.r.data(), tid_r, ds_mem, ds_r);
	// read particle velocities
	ds_v.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dset_v.read(sample.v.data(), tid_v, ds_mem, ds_v);
    }
    catch (H5::Exception const&) {
	throw exception("failed to read sample from HDF5 trajectory input file");
    }
}

/**
 * write phase space sample
 */
template <typename sample_type>
void trajectory::write(sample_type const& sample, double time)
{
    enum { dimension = sample_type::position_vector::static_size };
    typedef typename sample_type::position_vector::value_type position_value_type;
    typedef typename sample_type::velocity_vector::value_type velocity_value_type;

    H5::DataType const tid_r(H5xx::ctype<position_value_type>::type);
    H5::DataType const tid_v(H5xx::ctype<velocity_value_type>::type);
    H5::DataType const tid_t(H5::PredType::NATIVE_DOUBLE);
    unsigned int const npart = sample.r.size();

    H5::DataSet dset_t, dset_r, dset_v;
    try {
	H5::Group root(m_file.openGroup("trajectory"));
	dset_r = root.openDataSet("r");
	dset_v = root.openDataSet("v");
	dset_t = root.openDataSet("t");
    }
    catch (H5::Exception const&) {
	H5::Group root(m_file.createGroup("trajectory"));

	// vector sample file dataspace
	hsize_t dim_vector[3] = { 0, npart, dimension };
	hsize_t max_dim_vector[3] = { H5S_UNLIMITED, npart, dimension };
	H5::DataSpace ds_vector_file(3, dim_vector, max_dim_vector);
	// enable dataset chunking with GZIP compression
	H5::DSetCreatPropList cparms;
	hsize_t chunk_dim[3] = { 1, npart, dimension };
	cparms.setChunk(3, chunk_dim);
	cparms.setDeflate(6);

	dset_r = root.createDataSet("r", tid_r, ds_vector_file, cparms);
	dset_v = root.createDataSet("v", tid_v, ds_vector_file, cparms);

	// scalar sample file dataspace
	hsize_t dim_scalar[1] = { 0 };
	hsize_t max_dim_scalar[1] = { H5S_UNLIMITED };
	H5::DataSpace ds_scalar_file(1, dim_scalar, max_dim_scalar);
	hsize_t chunk_scalar[1] = { 1 };
	cparms.setChunk(1, chunk_scalar);

	dset_t = root.createDataSet("t", tid_t, ds_scalar_file, cparms);
    }

    // extend vector sample file dataspace
    H5::DataSpace ds_vector_file(dset_r.getSpace());
    hsize_t dim_vector[3];
    ds_vector_file.getSimpleExtentDims(dim_vector);
    hsize_t count_vector[3]  = { 1, 1, 1 };
    hsize_t start_vector[3]  = { dim_vector[0], 0, 0 };
    hsize_t stride_vector[3] = { 1, 1, 1 };
    hsize_t block_vector[3]  = { 1, npart, dimension };
    dim_vector[0]++;
    ds_vector_file.setExtentSimple(3, dim_vector);
    ds_vector_file.selectHyperslab(H5S_SELECT_SET, count_vector, start_vector, stride_vector, block_vector);

    // extend scalar sample file dataspace
    hsize_t count_scalar[1]  = { 1 };
    hsize_t start_scalar[1]  = { start_vector[0] };
    hsize_t stride_scalar[1] = { 1 };
    hsize_t block_scalar[1]  = { 1 };
    hsize_t dim_scalar[1] = { dim_vector[0] };
    H5::DataSpace ds_scalar_file(1, dim_scalar);
    ds_scalar_file.selectHyperslab(H5S_SELECT_SET, count_scalar, start_scalar, stride_scalar, block_scalar);

    // vector sample memory dataspace
    hsize_t dim_mem[2] = { npart, dimension };
    H5::DataSpace ds_vector_sample(2, dim_mem);

    // write periodically extended particle coordinates
    dset_r.extend(dim_vector);
    dset_r.write(sample.r.data(), tid_r, ds_vector_sample, ds_vector_file);
    // write particle velocities
    dset_v.extend(dim_vector);
    dset_v.write(sample.v.data(), tid_v, ds_vector_sample, ds_vector_file);
    // write simulation time
    dset_t.extend(dim_scalar);
    dset_t.write(&time, H5::PredType::NATIVE_DOUBLE, H5S_SCALAR, ds_scalar_file);
}

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_TRAJECTORY_HPP */
