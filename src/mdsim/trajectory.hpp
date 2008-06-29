/* MD simulation trajectory writer
 *
 * Copyright (C) 2008  Peter Colberg
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

#ifndef MDSIM_TRAJECTORY_HPP
#define MDSIM_TRAJECTORY_HPP

#include <H5Cpp.h>
#include <algorithm>
#include <assert.h>
#include <boost/array.hpp>
#include <cuda_wrapper.hpp>
#include <string>
#include "H5param.hpp"
#include "block.hpp"
#include "log.hpp"
#include "exception.hpp"


namespace mdsim {

template <unsigned dimension, typename T, typename U, bool writer = true>
class trajectory;

/**
 * trajectory file writer
 */
template <unsigned dimension, typename T, typename U>
class trajectory<dimension, T, U, true>
{
public:
    /** vector sample vector in page-locked host memory */
    typedef cuda::host::vector<U> vector_type;

public:
    trajectory(block_param<dimension, T> const& param);
    /** create HDF5 trajectory output file */
    void open(std::string const& filename, unsigned int const& npart);
    /** close HDF5 trajectory output file */
    void close();
    /** dump global simulation parameters to HDF5 file */
    trajectory<dimension, T, U, true>& operator<<(H5param const& param);
    /** write phase space sample */
    void sample(vector_type const& r, vector_type const& v, unsigned int const& npart, float const& timestep);

private:
    /** block algorithm parameters */
    block_param<dimension, T> param;
    /** number of samples */
    uint64_t samples_;

    /** HDF5 trajectory output file */
    H5::H5File file_;
    /** trajectory datasets for particle coordinates and velocities */
    boost::array<H5::DataSet, 3> dataset_;
    /** memory dataspace for a single coordinates or velocities sample */
    H5::DataSpace ds_mem_;
    /** file dataspace for a single coordinates or velocities sample */
    H5::DataSpace ds_file_;
    /** file dataspace for simulation time */
    H5::DataSpace ds_scalar_;
};

/**
 * trajectory file reader
 */
template <unsigned dimension, typename T, typename U>
class trajectory<dimension, T, U, false>
{
public:
    /** sample vector in page-locked host memory */
    typedef cuda::host::vector<U> vector_type;

public:
    trajectory();
    /** open HDF5 trajectory input file */
    void open(std::string const& filename);
    /** close HDF5 trajectory input file */
    void close();
    /** read global simulation parameters */
    void read(H5param& param);
    /** read phase space sample */
    void read(vector_type& r, vector_type& v, int64_t index);

private:
    /** HDF5 trajectory input file */
    H5::H5File file;
};

template <unsigned dimension, typename T, typename U>
trajectory<dimension, T, U, true>::trajectory(block_param<dimension, T> const& param) : param(param), samples_(0)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif
}

/**
 * create HDF5 trajectory output file
 */
template <unsigned dimension, typename T, typename U>
void trajectory<dimension, T, U, true>::open(std::string const& filename, unsigned int const& npart)
{
    // create trajectory output file
    LOG("write trajectories to file: " << filename);
    try {
	// truncate existing file
	file_ = H5::H5File(filename, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create trajectory output file");
    }

    // modify dataset creation properties to enable chunking
    H5::DSetCreatPropList cparms;
    hsize_t chunk_dim[3] = { 1, npart, dimension };
    cparms.setChunk(3, chunk_dim);

    // create datasets
    hsize_t dim[3] = { 0, npart, dimension };
    hsize_t max_dim[3] = { H5S_UNLIMITED, npart, dimension };
    ds_file_ = H5::DataSpace(3, dim, max_dim);
    H5::Group root(file_.createGroup("trajectory"));
    dataset_[0] = root.createDataSet("positions", H5::PredType::NATIVE_FLOAT, ds_file_, cparms);
    dataset_[1] = root.createDataSet("velocities", H5::PredType::NATIVE_FLOAT, ds_file_, cparms);

    hsize_t dim_mem[2] = { npart, sizeof(U) / sizeof(float) };
    ds_mem_ = H5::DataSpace(2, dim_mem);

    hsize_t count_mem[2]  = { npart, 1 };
    hsize_t start_mem[2]  = { 0, 0 };
    hsize_t stride_mem[2] = { 1, 1 };
    hsize_t block_mem[2]  = { 1, dimension };
    ds_mem_.selectHyperslab(H5S_SELECT_SET, count_mem, start_mem, stride_mem, block_mem);

    hsize_t chunk_scalar[1] = { 1 };
    cparms.setChunk(1, chunk_scalar);

    hsize_t dim_scalar[1] = { 0 };
    hsize_t max_dim_scalar[1] = { H5S_UNLIMITED };
    ds_scalar_ = H5::DataSpace(1, dim_scalar, max_dim_scalar);
    dataset_[2] = root.createDataSet("time", H5::PredType::NATIVE_FLOAT, ds_scalar_, cparms);
}

/**
 * close HDF5 trajectory output file
 */
template <unsigned dimension, typename T, typename U>
void trajectory<dimension, T, U, true>::close()
{
    try {
	file_.close();
    }
    catch (H5::Exception const& e) {
	throw exception("failed to close HDF5 trajectory input file");
    }
}

/**
 * dump global simulation parameters to HDF5 file
 */
template <unsigned dimension, typename T, typename U>
trajectory<dimension, T, U, true>& trajectory<dimension, T, U, true>::operator<<(H5param const& param)
{
    param.write(file_.createGroup("/parameters"));
    return *this;
}

/**
 * write phase space sample
 */
template <unsigned dimension, typename T, typename U>
void trajectory<dimension, T, U, true>::sample(vector_type const& r, vector_type const& v, unsigned int const& npart, float const& timestep)
{
    assert(r.size() == npart);
    assert(v.size() == npart);

    // absolute simulation time
    const float time_ = samples_ * timestep;

    hsize_t dim[3] = { samples_ + 1, npart, dimension };
    ds_file_.setExtentSimple(3, dim);

    hsize_t count[3]  = { 1, npart, 1 };
    hsize_t start[3]  = { samples_, 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, 1, dimension };

    ds_file_.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    hsize_t dim_scalar[1] = { samples_ + 1 };
    ds_scalar_.setExtentSimple(1, dim_scalar);

    hsize_t count_scalar[1]  = { 1 };
    hsize_t start_scalar[1]  = { samples_ };
    hsize_t stride_scalar[1] = { 1 };
    hsize_t block_scalar[1]  = { 1 };

    ds_scalar_.selectHyperslab(H5S_SELECT_SET, count_scalar, start_scalar, stride_scalar, block_scalar);

    // write periodically reduced particle coordinates
    dataset_[0].extend(dim);
    dataset_[0].write(r.data(), H5::PredType::NATIVE_FLOAT, ds_mem_, ds_file_);
    // write particle velocities
    dataset_[1].extend(dim);
    dataset_[1].write(v.data(), H5::PredType::NATIVE_FLOAT, ds_mem_, ds_file_);
    // write simulation time
    dataset_[2].extend(dim_scalar);
    dataset_[2].write(&time_, H5::PredType::NATIVE_FLOAT, H5S_SCALAR, ds_scalar_);

    samples_++;
}


template <unsigned dimension, typename T, typename U>
trajectory<dimension, T, U, false>::trajectory()
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif
}

/**
 * open HDF5 trajectory input file
 */
template <unsigned dimension, typename T, typename U>
void trajectory<dimension, T, U, false>::open(std::string const& filename)
{
    LOG("read trajectory file: " << filename);
    try {
	file = H5::H5File(filename, H5F_ACC_RDONLY);
    }
    catch (H5::Exception const& e) {
	throw exception("failed to open HDF5 trajectory input file");
    }
}

/**
 * close HDF5 trajectory input file
 */
template <unsigned dimension, typename T, typename U>
void trajectory<dimension, T, U, false>::close()
{
    try {
	file.close();
    }
    catch (H5::Exception const& e) {
	throw exception("failed to close HDF5 trajectory input file");
    }
}

/**
 * read global simulation parameters
 */
template <unsigned dimension, typename T, typename U>
void trajectory<dimension, T, U, false>::read(H5param& param)
{
    param.read(file.openGroup("/parameters"));
}

/**
 * read phase space sample
 */
template <unsigned dimension, typename T, typename U>
void trajectory<dimension, T, U, false>::read(vector_type& r, vector_type& v, int64_t index)
{
    try {
	// open phase space coordinates datasets
	H5::Group root(file.openGroup("trajectory"));
	H5::DataSet dataset_r(root.openDataSet("positions"));
	H5::DataSpace ds_r(dataset_r.getSpace());
	H5::DataSet dataset_v(root.openDataSet("velocities"));
	H5::DataSpace ds_v(dataset_v.getSpace());

	// validate dataspace extents
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

	// retrieve dataspace dimensions
	hsize_t dim_r[3];
	ds_r.getSimpleExtentDims(dim_r);
	hsize_t dim_v[3];
	ds_v.getSimpleExtentDims(dim_v);

	if (!std::equal(dim_r, dim_r + 3, dim_v)) {
	    throw exception("trajectory and velocity dataspace dimensions differ");
	}

	// number of samples
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

	assert(r.size() == npart);
	assert(v.size() == npart);

	// read sample from dataset
	hsize_t dim_mem[2] = { npart, sizeof(U) / sizeof(float) };
	H5::DataSpace ds_mem(2, dim_mem);

	hsize_t count_mem[2]  = { npart, 1 };
	hsize_t start_mem[2]  = { 0, 0 };
	hsize_t stride_mem[2] = { 1, 1 };
	hsize_t block_mem[2]  = { 1, dimension };
	ds_mem.selectHyperslab(H5S_SELECT_SET, count_mem, start_mem, stride_mem, block_mem);

	hsize_t count[3]  = { 1, npart, 1 };
	hsize_t start[3]  = { index, 0, 0 };
	hsize_t stride[3] = { 1, 1, 1 };
	hsize_t block[3]  = { 1, 1, dimension };

	// coordinates
	ds_r.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dataset_r.read(r.data(), H5::PredType::NATIVE_FLOAT, ds_mem, ds_r);
	// velocities
	ds_v.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dataset_v.read(v.data(), H5::PredType::NATIVE_FLOAT, ds_mem, ds_v);
    }
    catch (H5::Exception const& e) {
	throw exception("failed to read from HDF5 trajectory input file");
    }
}

} // namespace mdsim

#endif /* ! MDSIM_TRAJECTORY_HPP */
