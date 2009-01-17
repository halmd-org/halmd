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
    void read(sample_type& sample, ssize_t index);
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
    void read_vector_sample(H5::DataSet dset, sample_type& sample, ssize_t& index);
    template <typename sample_type>
    void read_scalar_sample(H5::DataSet dset, sample_type& sample, ssize_t& index);
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
void trajectory::read(sample_type& sample, ssize_t index)
{
    std::vector<H5::DataSet> dset_r;
    std::vector<H5::DataSet> dset_v;
    H5::DataSet dset_t;

    try {
	H5XX_NO_AUTO_PRINT(H5::Exception);
	H5::Group root(m_file.openGroup("trajectory"));
	try {
	    // binary mixture
	    H5::Group a(root.openGroup("A"));
	    dset_r.push_back(a.openDataSet("r"));
	    dset_v.push_back(a.openDataSet("v"));
	    H5::Group b(root.openGroup("B"));
	    dset_r.push_back(b.openDataSet("r"));
	    dset_v.push_back(b.openDataSet("v"));
	}
	catch (H5::GroupIException const&) {
	    try {
		// backwards compatibility with r:R:v:t format
		// 	 r = reduced single- or double-precision positions,
		//   R = extended single- or double-precision positions,
		//   v = single- or double-precision velocities
		//   t = single- or double-precision simulation time
		dset_r.push_back(root.openDataSet("R"));
		LOG_WARNING("detected obsolete trajectory file format");
	    }
	    catch (H5::GroupIException const&)
	    {
		// new-style r:v:t format
		//   r = extended double-precision positions,
		//   v = single- or double-precision velocities
		//   t = double-precision simulation time
		dset_r.push_back(root.openDataSet("r"));
	    }
	    // backwards compatibility with r:R:v:t format
	    if (dset_r.back().getDataType() == H5::PredType::NATIVE_FLOAT) {
		// use reduced positions if extended positions are single-precision
		dset_r.pop_back();
		dset_r.push_back(root.openDataSet("r"));
		LOG_WARNING("falling back to reduced particle position sample");
	    }
	    dset_v.push_back(root.openDataSet("v"));
	}
	dset_t = root.openDataSet("t");
    }
    catch (H5::Exception const&) {
	throw exception("failed to open HDF5 trajectory datasets");
    }

    for (size_t i = 0; i < dset_r.size(); ++i) {
	read_vector_sample(dset_r[i], sample[i].r, index);
	read_vector_sample(dset_v[i], sample[i].v, index);
    }
    read_scalar_sample(dset_t, sample.time, index);

    LOG("resuming from trajectory sample at offset " << index << " with time " << sample.time);

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
void trajectory::read_vector_sample(H5::DataSet dset, sample_type& sample, ssize_t& index)
{
    typedef typename sample_type::value_type::value_type float_type;
    enum { dimension = sample_type::value_type::static_size };

    H5::DataType const tid(H5xx::ctype<float_type>::type);
    H5::DataSpace ds(dset.getSpace());

    if (!ds.isSimple()) {
	throw exception("HDF5 vector dataspace is not a simple dataspace");
    }
    if (ds.getSimpleExtentNdims() != 3) {
	throw exception("HDF5 vector dataspace has invalid dimensionality");
    }

    hsize_t dim[3];
    ds.getSimpleExtentDims(dim);

    ssize_t const len = dim[0];
    if ((index >= len) || ((-index) > len)) {
	throw exception("trajectory input sample number out of bounds");
    }
    index = (index < 0) ? (index + len) : index;

    size_t const size = dim[1];
    if (size < 1) {
	throw exception("trajectory input file has invalid number of particles");
    }
    if (dim[2] != dimension) {
	throw exception("trajectory input file has invalid coordinate dimension");
    }

    try {
	sample.resize(size);
    }
    catch (std::bad_alloc const&) {
	throw exception("failed to allocate memory for phase space sample vector");
    }

    hsize_t dim_sample[2] = { size, dimension };
    H5::DataSpace ds_sample(2, dim_sample);

    hsize_t count[3]  = { 1, size, 1 };
    hsize_t start[3]  = { index, 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, 1, dimension };
    ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    try {
	H5XX_NO_AUTO_PRINT(H5::Exception);
	dset.read(sample.data(), tid, ds_sample, ds);
    }
    catch (H5::Exception const&) {
	throw exception("failed to read vector sample from HDF5 trajectory input file");
    }
}

template <typename sample_type>
void trajectory::read_scalar_sample(H5::DataSet dset, sample_type& sample, ssize_t& index)
{
    H5::DataType const tid(H5xx::ctype<sample_type>::type);
    H5::DataSpace ds(dset.getSpace());

    if (!ds.isSimple()) {
	throw exception("HDF5 scalar dataspace is not a simple dataspace");
    }
    if (ds.getSimpleExtentNdims() != 1) {
	throw exception("HDF5 scalar dataspace has invalid dimensionality");
    }

    hsize_t dim[1];
    ds.getSimpleExtentDims(dim);

    ssize_t const len = dim[0];
    if ((index >= len) || ((-index) > len)) {
	throw exception("trajectory input sample number out of bounds");
    }
    index = (index < 0) ? (index + len) : index;

    hsize_t count[1]  = { 1 };
    hsize_t start[1]  = { index };
    hsize_t stride[1] = { 1 };
    hsize_t block[1]  = { 1 };
    ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    try {
	H5XX_NO_AUTO_PRINT(H5::Exception);
	dset.read(&sample, tid, H5S_SCALAR, ds);
    }
    catch (H5::Exception const&) {
	throw exception("failed to read scalar sample from HDF5 trajectory input file");
    }
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
