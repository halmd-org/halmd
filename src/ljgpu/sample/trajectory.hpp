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
#include <ljgpu/sample/sample.hpp>
#include <ljgpu/util/H5param.hpp>
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
    /** read phase space sample */
    template <typename sample_vector>
    void read(sample_vector& r, sample_vector& v, int64_t index);
    /** write phase space sample */
    template <typename sample_type>
    void write(sample_type const& sample, double time);
    /** flush HDF5 output file to disk */
    void flush();
    /** returns HDF5 parameter group */
    H5param attrs();

private:
    /** create HDF5 datasets */
    template <int dimension>
    void create_datasets(trajectory_gpu_sample<vector<float, dimension> > const&,
			 H5::DSetCreatPropList const& cparms);
    template <int dimension>
    void create_datasets(trajectory_host_sample<vector<double, dimension> > const&,
			 H5::DSetCreatPropList const& cparms);
    /** write HDF5 datasets */
    template <int dimension>
    void write_datasets(trajectory_gpu_sample<vector<float, dimension> > const& sample,
			hsize_t const* dim);
    template <int dimension>
    void write_datasets(trajectory_host_sample<vector<double, dimension> > const& sample,
			hsize_t const* dim);

private:
    /** HDF5 trajectory output file */
    H5::H5File m_file;
    /** trajectory datasets for particle coordinates and velocities */
    boost::array<H5::DataSet, 4> m_dataset;
    /** memory dataspace for a single coordinates or velocities sample */
    H5::DataSpace m_ds_mem;
    /** file dataspace for a single coordinates or velocities sample */
    H5::DataSpace m_ds_file;
    /** file dataspace for simulation time */
    H5::DataSpace m_ds_scalar;
    /** floating-point data type */
    H5::DataType m_tid;
    /** sample index */
    int64_t m_index;
};

/**
 * read phase space sample
 */
template <typename sample_vector>
void trajectory::read(sample_vector& r, sample_vector& v, int64_t index)
{
    enum { dimension = sample_vector::value_type::static_size };
    typedef typename sample_vector::value_type::value_type float_type;

    H5::DataSet dataset_r, dataset_v;
    try {
	H5::Group root(m_file.openGroup("trajectory"));
	try {
	    // read reduced particle coordinates, if available
	    dataset_r = root.openDataSet("r");
	}
	catch (H5::Exception const&) {
	    LOG_WARNING("falling back to extended particle coordinates");
	    dataset_r = root.openDataSet("R");
	}
	dataset_v = root.openDataSet("v");
    }
    catch (H5::Exception const&) {
	throw exception("failed to open HDF5 trajectory datasets");
    }

    // validate dataspace extents
    H5::DataSpace ds_r(dataset_r.getSpace());
    H5::DataSpace ds_v(dataset_v.getSpace());
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
	r.resize(npart);
	v.resize(npart);
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

    if (boost::is_same<float_type, double>::value)
	m_tid = H5::PredType::NATIVE_DOUBLE;
    else
	m_tid = H5::PredType::NATIVE_FLOAT;

    try {
	// read periodically reduced particle positions
	ds_r.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dataset_r.read(r.data(), m_tid, ds_mem, ds_r);
	// read particle velocities
	ds_v.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dataset_v.read(v.data(), m_tid, ds_mem, ds_v);
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
    enum { dimension = sample_type::sample_vector::value_type::static_size };
    typedef typename sample_type::sample_vector::value_type::value_type float_type;

    // create datasets
    if (!m_index) {
	// create position and velocity datasets
	unsigned int npart = sample.R.size();
	hsize_t dim[3] = { 0, npart, dimension };
	hsize_t max_dim[3] = { H5S_UNLIMITED, npart, dimension };
	m_ds_file = H5::DataSpace(3, dim, max_dim);
	// modify dataset creation properties to enable chunking
	H5::DSetCreatPropList cparms;
	hsize_t chunk_dim[3] = { 1, npart, dimension };
	cparms.setChunk(3, chunk_dim);
	// use GZIP compression
	cparms.setDeflate(6);
	create_datasets<dimension>(sample, cparms);

	// create time dataset
	H5::Group root(m_file.openGroup("trajectory"));
	hsize_t dim_mem[2] = { npart, dimension };
	m_ds_mem = H5::DataSpace(2, dim_mem);
	hsize_t chunk_scalar[1] = { 1 };
	cparms.setChunk(1, chunk_scalar);
	hsize_t dim_scalar[1] = { 0 };
	hsize_t max_dim_scalar[1] = { H5S_UNLIMITED };
	m_ds_scalar = H5::DataSpace(1, dim_scalar, max_dim_scalar);
	m_dataset[0] = root.createDataSet("t", H5::PredType::NATIVE_DOUBLE, m_ds_scalar, cparms);
    }

    hsize_t dim[3];
    m_ds_file.getSimpleExtentDims(dim);

    hsize_t count[3]  = { 1, 1, 1 };
    hsize_t start[3]  = { dim[0], 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, dim[1], dim[2] };

    dim[0]++;
    m_ds_file.setExtentSimple(3, dim);
    m_ds_file.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    write_datasets<dimension>(sample, dim);

    hsize_t dim_scalar[1];
    m_ds_scalar.getSimpleExtentDims(dim_scalar);

    hsize_t count_scalar[1]  = { 1 };
    hsize_t start_scalar[1]  = { dim_scalar[0] };
    hsize_t stride_scalar[1] = { 1 };
    hsize_t block_scalar[1]  = { 1 };

    dim_scalar[0]++;
    m_ds_scalar.setExtentSimple(1, dim_scalar);
    m_ds_scalar.selectHyperslab(H5S_SELECT_SET, count_scalar, start_scalar, stride_scalar, block_scalar);

    // write simulation time
    m_dataset[0].extend(dim_scalar);
    m_dataset[0].write(&time, H5::PredType::NATIVE_DOUBLE, H5S_SCALAR, m_ds_scalar);

    m_index++;
}

template <int dimension>
void trajectory::create_datasets(trajectory_gpu_sample<vector<float, dimension> > const&,
				 H5::DSetCreatPropList const& cparms)
{
    m_tid = H5::PredType::NATIVE_FLOAT;
    try {
	H5::Group root(m_file.createGroup("trajectory"));
	m_dataset[1] = root.createDataSet("r", m_tid, m_ds_file, cparms);
	m_dataset[2] = root.createDataSet("R", m_tid, m_ds_file, cparms);
	m_dataset[3] = root.createDataSet("v", m_tid, m_ds_file, cparms);
    }
    catch (H5::Exception const&) {
	throw exception("failed to create HDF5 trajectory datasets");
    }
}

template <int dimension>
void trajectory::create_datasets(trajectory_host_sample<vector<double, dimension> > const&,
				 H5::DSetCreatPropList const& cparms)
{
    m_tid = H5::PredType::NATIVE_DOUBLE;
    try {
	H5::Group root(m_file.createGroup("trajectory"));
	m_dataset[1] = root.createDataSet("R", m_tid, m_ds_file, cparms);
	m_dataset[2] = root.createDataSet("v", m_tid, m_ds_file, cparms);
    }
    catch (H5::Exception const&) {
	throw exception("failed to create HDF5 trajectory datasets");
    }
}

template <int dimension>
void trajectory::write_datasets(trajectory_gpu_sample<vector<float, dimension> > const& sample,
				hsize_t const* dim)
{
    // write periodically reduced particle coordinates
    m_dataset[1].extend(dim);
    m_dataset[1].write(sample.r.data(), m_tid, m_ds_mem, m_ds_file);
    // write periodically extended particle coordinates
    m_dataset[2].extend(dim);
    m_dataset[2].write(sample.R.data(), m_tid, m_ds_mem, m_ds_file);
    // write particle velocities
    m_dataset[3].extend(dim);
    m_dataset[3].write(sample.v.data(), m_tid, m_ds_mem, m_ds_file);
}

template <int dimension>
void trajectory::write_datasets(trajectory_host_sample<vector<double, dimension> > const& sample,
				hsize_t const* dim)
{
    // write periodically extended particle coordinates
    m_dataset[1].extend(dim);
    m_dataset[1].write(sample.R.data(), m_tid, m_ds_mem, m_ds_file);
    // write particle velocities
    m_dataset[2].extend(dim);
    m_dataset[2].write(sample.v.data(), m_tid, m_ds_mem, m_ds_file);
}

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_TRAJECTORY_HPP */
