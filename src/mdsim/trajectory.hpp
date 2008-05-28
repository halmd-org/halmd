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
#include <string>
#include "H5param.hpp"
#include "block.hpp"
#include "log.hpp"
#include "exception.hpp"


namespace mdsim {

/**
 * trajectory file reader or writer
 */
template <unsigned dimension, typename T, bool writer = true>
class trajectory;


/**
 * trajectory file writer
 */
template <unsigned dimension, typename T>
class trajectory<dimension, T, true>
{
public:
    trajectory(block_param<dimension, T> const& param);
    /** create HDF5 trajectory output file */
    void open(std::string const& filename, unsigned int const& npart);
    /** close HDF5 trajectory output file */
    void close();
    /** dump global simulation parameters to HDF5 file */
    trajectory<dimension, T, true>& operator<<(H5param const& param);
    /** write phase space sample */
    void sample(std::vector<T> const& r, std::vector<T> const& v, unsigned int const& npart);

private:
    /** block algorithm parameters */
    block_param<dimension, T> param;
    /** number of samples */
    uint64_t samples_;

    /** HDF5 trajectory output file */
    H5::H5File file_;
    /** trajectory datasets for particle coordinates and velocities */
    boost::array<H5::DataSet, 2> dataset_;
    /** memory dataspace for a single coordinates or velocities sample */
    H5::DataSpace ds_mem_;
    /** file dataspace for a single coordinates or velocities sample */
    H5::DataSpace ds_file_;
};


template <unsigned dimension, typename T>
trajectory<dimension, T, true>::trajectory(block_param<dimension, T> const& param) : param(param), samples_(0)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif
}

/**
 * create HDF5 trajectory output file
 */
template <unsigned dimension, typename T>
void trajectory<dimension, T, true>::open(std::string const& filename, unsigned int const& npart)
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

    // create datasets
    hsize_t dim[3] = { param.max_samples(), npart, dimension };
    ds_file_ = H5::DataSpace(3, dim);
    dataset_[0] = file_.createDataSet("positions", H5::PredType::NATIVE_DOUBLE, ds_file_);
    dataset_[1] = file_.createDataSet("velocities", H5::PredType::NATIVE_DOUBLE, ds_file_);

    hsize_t dim_mem[2] = { npart, dimension };
    ds_mem_ = H5::DataSpace(2, dim_mem);
}

/**
 * close HDF5 trajectory output file
 */
template <unsigned dimension, typename T>
void trajectory<dimension, T, true>::close()
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
template <unsigned dimension, typename T>
trajectory<dimension, T, true>& trajectory<dimension, T, true>::operator<<(H5param const& param)
{
    param.write(file_.createGroup("/parameters"));
    return *this;
}

/**
 * write phase space sample
 */
template <unsigned dimension, typename T>
void trajectory<dimension, T, true>::sample(std::vector<T> const& r, std::vector<T> const& v, unsigned int const& npart)
{
    if (samples_ >= param.max_samples())
	return;

    assert(r.size() == npart);
    assert(v.size() == npart);

    hsize_t count[3]  = { 1, npart, 1 };
    hsize_t start[3]  = { samples_, 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, 1, dimension };

    ds_file_.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    // write particle positions
    dataset_[0].write(r.data(), H5::PredType::NATIVE_DOUBLE, ds_mem_, ds_file_);
    // write particle velocities
    dataset_[1].write(v.data(), H5::PredType::NATIVE_DOUBLE, ds_mem_, ds_file_);

    samples_++;
}


/**
 * trajectory file reader
 */
template <unsigned dimension, typename T>
class trajectory<dimension, T, false>
{
public:
    trajectory();
    /** open HDF5 trajectory input file */
    void open(std::string const& filename);
    /** close HDF5 trajectory input file */
    void close();
    /** read global simulation parameters */
    void read(H5param& param);
    /** read phase space sample */
    void read(std::vector<T>& r, std::vector<T>& v, int64_t index);

private:
    /** HDF5 trajectory input file */
    H5::H5File file;
};

template <unsigned dimension, typename T>
trajectory<dimension, T, false>::trajectory()
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif
}

/**
 * open HDF5 trajectory input file
 */
template <unsigned dimension, typename T>
void trajectory<dimension, T, false>::open(std::string const& filename)
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
template <unsigned dimension, typename T>
void trajectory<dimension, T, false>::close()
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
template <unsigned dimension, typename T>
void trajectory<dimension, T, false>::read(H5param& param)
{
    param.read(file.openGroup("/parameters"));
}

/**
 * read phase space sample
 */
template <unsigned dimension, typename T>
void trajectory<dimension, T, false>::read(std::vector<T>& r, std::vector<T>& v, int64_t index)
{
    try {
	// open phase space coordinates datasets
	H5::DataSet dataset_r(file.openDataSet("positions"));
	H5::DataSpace ds_r(dataset_r.getSpace());
	H5::DataSet dataset_v(file.openDataSet("velocities"));
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
	hsize_t dim_mem[2] = { npart, dimension };
	H5::DataSpace ds_mem(2, dim_mem);

	hsize_t count[3]  = { 1, npart, 1 };
	hsize_t start[3]  = { index, 0, 0 };
	hsize_t stride[3] = { 1, 1, 1 };
	hsize_t block[3]  = { 1, 1, dimension };

	// read particle positions
	ds_r.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dataset_r.read(r.data(), H5::PredType::NATIVE_DOUBLE, ds_mem, ds_r);
	// read particle velocities
	ds_v.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dataset_v.read(v.data(), H5::PredType::NATIVE_DOUBLE, ds_mem, ds_v);
    }
    catch (H5::Exception const& e) {
	throw exception("failed to read from HDF5 trajectory input file");
    }
}

} // namespace mdsim

#endif /* ! MDSIM_TRAJECTORY_HPP */
