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
#include <ljgpu/sample/H5param.hpp>
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
    operator H5param() { return m_file; }

private:
    template <typename sample_type>
    H5::DataSet create_vector_dataset(H5::Group node, char const* name, sample_type const& sample);
    template <typename sample_type>
    H5::DataSet create_scalar_dataset(H5::Group node, char const* name, sample_type const&);
    template <typename sample_type>
    void write_vector_sample(H5::DataSet dset, sample_type const& sample);
    template <typename sample_type>
    void write_scalar_sample(H5::DataSet dset, sample_type const& sample);

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
	H5XX_NO_AUTO_PRINT(H5::FileIException);
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
    catch (H5::FileIException const&) {
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
	sample[0 /* FIXME */].r.resize(npart);
	sample[0 /* FIXME */].v.resize(npart);
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
	H5XX_NO_AUTO_PRINT(H5::Exception);
	// read periodically reduced particle positions
	ds_r.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dset_r.read(sample[0 /* FIXME */].r.data(), tid_r, ds_mem, ds_r);
	// read particle velocities
	ds_v.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dset_v.read(sample[0 /* FIXME */].v.data(), tid_v, ds_mem, ds_v);
    }
    catch (H5::Exception const&) {
	throw exception("failed to read sample from HDF5 trajectory input file");
    }

    try {
	H5XX_NO_AUTO_PRINT(H5::Exception);
	sample.box = H5param(*this)["mdsim"]["box_length"].as<float>();
    }
    catch (H5::Exception const&) {
	throw exception("failed to read simulation box length from HDF5 file");
    }
}

/**
 * write phase space sample
 */
template <typename sample_type>
void trajectory::write(sample_type const& sample, double time)
{
    boost::array<H5::DataSet, 2> dset_r;
    boost::array<H5::DataSet, 2> dset_v;
    H5::DataSet dset_t;

    try {
	H5XX_NO_AUTO_PRINT(H5::FileIException);
	H5::Group root(m_file.openGroup("trajectory"));

	if (sample[1].r.size()) {
	    // binary mixture
	    H5::Group a(root.openGroup("A"));
	    dset_r[0] = a.openDataSet("r");
	    dset_v[0] = a.openDataSet("v");
	    H5::Group b(root.openGroup("B"));
	    dset_r[1] = b.openDataSet("r");
	    dset_v[1] = b.openDataSet("v");
	}
	else {
	    dset_r[0] = root.openDataSet("r");
	    dset_v[0] = root.openDataSet("v");
	}
	dset_t = root.openDataSet("t");
    }
    catch (H5::FileIException const&) {
	H5::Group root(m_file.createGroup("trajectory"));

	if (sample[1].r.size()) {
	    // binary mixture
	    H5::Group a(root.createGroup("A"));
	    dset_r[0] = create_vector_dataset(a, "r", sample[0].r);
	    dset_v[0] = create_vector_dataset(a, "v", sample[0].v);
	    H5::Group b(root.createGroup("B"));
	    dset_r[1] = create_vector_dataset(b, "r", sample[1].r);
	    dset_v[1] = create_vector_dataset(b, "v", sample[1].v);
	}
	else {
	    dset_r[0] = create_vector_dataset(root, "r", sample[0].r);
	    dset_v[0] = create_vector_dataset(root, "v", sample[0].v);
	}
	dset_t = create_scalar_dataset(root, "t", time);
    }

    write_vector_sample(dset_r[0], sample[0].r);
    write_vector_sample(dset_v[0], sample[0].v);

    if (sample[1].r.size()) {
	write_vector_sample(dset_r[1], sample[1].r);
	write_vector_sample(dset_v[1], sample[1].v);
    }

    write_scalar_sample(dset_t, time);
}

template <typename sample_type>
H5::DataSet trajectory::create_vector_dataset(H5::Group node, char const* name, sample_type const& sample)
{
    typedef typename sample_type::value_type::value_type float_type;
    enum { dimension = sample_type::value_type::static_size };

    H5::DataType const tid(H5xx::ctype<float_type>::type);
    size_t const size = sample.size();

    // vector sample file dataspace
    hsize_t dim[3] = { 0, size, dimension };
    hsize_t max_dim[3] = { H5S_UNLIMITED, size, dimension };
    H5::DataSpace ds(3, dim, max_dim);

    H5::DSetCreatPropList cparms;
    hsize_t chunk_dim[3] = { 1, size, dimension };
    cparms.setChunk(3, chunk_dim);
    // enable GZIP compression
    cparms.setDeflate(6);

    return node.createDataSet(name, tid, ds, cparms);
}

template <typename sample_type>
H5::DataSet trajectory::create_scalar_dataset(H5::Group node, char const* name, sample_type const&)
{
    H5::DataType const tid(H5xx::ctype<sample_type>::type);

    // scalar sample file dataspace
    hsize_t dim[1] = { 0 };
    hsize_t max_dim[1] = { H5S_UNLIMITED };
    H5::DataSpace ds(1, dim, max_dim);

    H5::DSetCreatPropList cparms;
    hsize_t chunk_dim[1] = { 1 };
    cparms.setChunk(1, chunk_dim);

    return node.createDataSet(name, tid, ds, cparms);
}

template <typename sample_type>
void trajectory::write_vector_sample(H5::DataSet dset, sample_type const& sample)
{
    typedef typename sample_type::value_type::value_type float_type;
    enum { dimension = sample_type::value_type::static_size };

    H5::DataType const tid(H5xx::ctype<float_type>::type);
    size_t const size = sample.size();

    // extend vector sample file dataspace
    H5::DataSpace ds(dset.getSpace());
    hsize_t dim[3];
    ds.getSimpleExtentDims(dim);
    hsize_t count[3]  = { 1, 1, 1 };
    hsize_t start[3]  = { dim[0], 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, size, dimension };
    dim[0]++;
    ds.setExtentSimple(3, dim);
    ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    // vector sample memory dataspace
    hsize_t dim_sample[2] = { size, dimension };
    H5::DataSpace ds_sample(2, dim_sample);

    dset.extend(dim);
    dset.write(sample.data(), tid, ds_sample, ds);
}

template <typename sample_type>
void trajectory::write_scalar_sample(H5::DataSet dset, sample_type const& sample)
{
    H5::DataType const tid(H5xx::ctype<sample_type>::type);

    // extend scalar sample file dataspace
    H5::DataSpace ds(dset.getSpace());
    hsize_t dim[1];
    ds.getSimpleExtentDims(dim);
    hsize_t count[1]  = { 1 };
    hsize_t start[1]  = { dim[0] };
    hsize_t stride[1] = { 1 };
    hsize_t block[1]  = { 1 };
    dim[0]++;
    ds.setExtentSimple(1, dim);
    ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    dset.extend(dim);
    dset.write(&sample, tid, H5S_SCALAR, ds);
}

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_TRAJECTORY_HPP */
