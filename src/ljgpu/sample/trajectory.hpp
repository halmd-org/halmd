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
#include <boost/noncopyable.hpp>
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
class trajectory : boost::noncopyable
{
public:
    /** io flags */
    enum openmode {
	in = 0x1,
	out = 0x2,
    };

public:
    trajectory() : m_is_open(false) {}
    /** create HDF5 trajectory output file */
    void open(std::string const& filename, openmode mode = in);
    /** close HDF5 trajectory output file */
    void close();
    /** returns true iff associated with HDF5 file */
    bool is_open() const { return m_is_open; }
    /** flush HDF5 output file to disk */
    void flush();

    /** read phase space sample */
    template <typename host_sample_type>
    void read(host_sample_type& sample, ssize_t index);
    /** write phase space sample */
    template <typename host_sample_type>
    void write(host_sample_type const& sample, double time);

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
    /** HDF5 data sets */
    std::vector<H5::DataSet> m_dset_r;
    std::vector<H5::DataSet> m_dset_v;
    H5::DataSet m_dset_t;
    /** true iff associated with HDF5 file */
    bool m_is_open;
};

/**
 * read phase space sample
 */
template <typename host_sample_type>
void trajectory::read(host_sample_type& sample, ssize_t index)
{
    typedef typename host_sample_type::value_type sample_type;
    typedef typename sample_type::position_sample_vector position_sample_vector;
    typedef typename sample_type::velocity_sample_vector velocity_sample_vector;
    typedef typename sample_type::position_sample_ptr position_sample_ptr;
    typedef typename sample_type::velocity_sample_ptr velocity_sample_ptr;

    std::vector<H5::DataSet> dset_r;
    std::vector<H5::DataSet> dset_v;
    H5::DataSet dset_t;
    std::vector<unsigned int> mpart;

    try {
	H5XX_NO_AUTO_PRINT(H5::Exception);
	H5::Group root(m_file.openGroup("trajectory"));

	// read particle numbers in binary mixture
	mpart = H5param(*this)["mdsim"]["particles"].as<std::vector<unsigned int> >();

	for (size_t i = 0; i < mpart.size(); ++i) {
	    H5::Group node(root);
	    if (mpart.size() > 1) {
		std::string name;
		name.push_back('A' + i);
		node = H5::Group(root.openGroup(name));
	    }
	    try {
		// backwards compatibility with r:R:v:t format
		//   r = reduced single- or double-precision positions,
		//   R = extended single- or double-precision positions,
		//   v = single- or double-precision velocities
		//   t = single- or double-precision simulation time
		dset_r.push_back(node.openDataSet("R"));
		LOG_WARNING("detected obsolete trajectory file format");
	    }
	    catch (H5::GroupIException const&)
	    {
		// new-style r:v:t format
		//   r = extended double-precision positions,
		//   v = single- or double-precision velocities
		//   t = double-precision simulation time
		dset_r.push_back(node.openDataSet("r"));
	    }
	    // backwards compatibility with r:R:v:t format
	    if (dset_r.back().getDataType() == H5::PredType::NATIVE_FLOAT) {
		// use reduced positions if extended positions are single-precision
		dset_r.pop_back();
		dset_r.push_back(node.openDataSet("r"));
		LOG_WARNING("falling back to reduced particle position sample");
	    }
	    dset_v.push_back(node.openDataSet("v"));
	}
	dset_t = root.openDataSet("t");
    }
    catch (H5::Exception const&) {
	throw exception("failed to open HDF5 trajectory datasets");
    }

    for (size_t i = 0; i < mpart.size(); ++i) {
	position_sample_ptr r(new position_sample_vector(mpart[i]));
	velocity_sample_ptr v(new velocity_sample_vector(mpart[i]));
	read_vector_sample(dset_r[i], *r, index);
	read_vector_sample(dset_v[i], *v, index);
	sample.push_back(sample_type(r, v));
    }
    double time;
    read_scalar_sample(dset_t, time, index);
    LOG("resuming from trajectory sample at offset " << index << " with time " << time);
}

/**
 * write phase space sample
 */
template <typename host_sample_type>
void trajectory::write(host_sample_type const& sample, double time)
{
    if (m_dset_r.empty()) {
	H5::Group root(m_file.createGroup("trajectory"));
	m_dset_r.reserve(sample.size());
	m_dset_v.reserve(sample.size());

	for (size_t i = 0; i < sample.size(); ++i) {
	    H5::Group node(root);
	    if (sample.size() > 1) {
		std::string name;
		name.push_back('A' + i);
		node = H5::Group(root.createGroup(name));
	    }
	    m_dset_r.push_back(create_vector_dataset(node, "r", *sample[i].r));
	    m_dset_v.push_back(create_vector_dataset(node, "v", *sample[i].v));
	}
	m_dset_t = create_scalar_dataset(root, "t", time);
    }

    assert(m_dset_r.size() == sample.size());
    assert(m_dset_v.size() == sample.size());

    for (size_t i = 0; i < sample.size(); ++i) {
	write_vector_sample(m_dset_r[i], *sample[i].r);
	write_vector_sample(m_dset_v[i], *sample[i].v);
    }
    write_scalar_sample(m_dset_t, time);
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
