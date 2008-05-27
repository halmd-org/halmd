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
    trajectory(H5param const& param);
    /** create HDF5 trajectory output file */
    void open(std::string const& filename);
    /** close HDF5 trajectory output file */
    void close();
    /** write global simulation parameters */
    trajectory<dimension, T, true>& operator<<(H5param const& param);
    /** write phase space sample */
    void sample(cuda::host::vector<T> const& r, cuda::host::vector<T> const& v);

private:
    /** HDF5 output file */
    H5::H5File file_;
    uint64_t npart_;
    uint64_t max_samples_;
    uint64_t samples_;
    boost::array<H5::DataSet, 2> dset_;
    H5::DataSpace ds_mem_;
    H5::DataSpace ds_file_;
};


template <unsigned dimension, typename T>
trajectory<dimension, T, true>::trajectory(H5param const& param) : npart_(param.particles()), max_samples_(param.max_samples()), samples_(0)
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
void trajectory<dimension, T, true>::open(std::string const& filename)
{
    // create trajectory output file
    LOG("trajectory output file: " << filename);
    try {
	// truncate existing file
	file_ = H5::H5File(filename, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create trajectory output file");
    }

    // create datasets
    hsize_t dim[3] = { max_samples_, npart_, dimension };
    ds_file_ = H5::DataSpace(3, dim);
    dset_[0] = file_.createDataSet("trajectory", H5::PredType::NATIVE_FLOAT, ds_file_);
    dset_[1] = file_.createDataSet("velocity", H5::PredType::NATIVE_FLOAT, ds_file_);

    hsize_t dim_mem[2] = { npart_, dimension };
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
 * write global simulation parameters
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
void trajectory<dimension, T, true>::sample(cuda::host::vector<T> const& r, cuda::host::vector<T> const& v)
{
    if (samples_ >= max_samples_)
	return;

    assert(r.size() == npart_);
    assert(v.size() == npart_);

    hsize_t count[3]  = { 1, npart_, 1 };
    hsize_t start[3]  = { samples_, 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, 1, dimension };

    ds_file_.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    // write periodically reduced particle coordinates
    dset_[0].write(r.data(), H5::PredType::NATIVE_FLOAT, ds_mem_, ds_file_);
    // write particle velocities
    dset_[1].write(v.data(), H5::PredType::NATIVE_FLOAT, ds_mem_, ds_file_);

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
    void read(cuda::host::vector<T>& r, cuda::host::vector<T>& v, int64_t index);

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
    LOG("resuming from trajectory file: " << filename);
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
void trajectory<dimension, T, false>::read(cuda::host::vector<T>& r, cuda::host::vector<T>& v, int64_t index)
{
    try {
	// open phase space coordinates datasets
	H5::DataSet dset_r(file.openDataSet("trajectory"));
	H5::DataSpace ds_r(dset_r.getSpace());
	H5::DataSet dset_v(file.openDataSet("velocity"));
	H5::DataSpace ds_v(dset_v.getSpace());

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

	LOG("reading trajectory sample at offset: " << index);

	// allocate memory for sample
	r.resize(npart);
	v.resize(npart);

	// read sample from dataset
	hsize_t dim_mem[2] = { npart, dimension };
	H5::DataSpace ds_mem(2, dim_mem);

	hsize_t count[3]  = { 1, npart, 1 };
	hsize_t start[3]  = { index, 0, 0 };
	hsize_t stride[3] = { 1, 1, 1 };
	hsize_t block[3]  = { 1, 1, dimension };

	// coordinates
	ds_r.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dset_r.read(r.data(), H5::PredType::NATIVE_FLOAT, ds_mem, ds_r);
	// velocities
	ds_v.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dset_v.read(v.data(), H5::PredType::NATIVE_FLOAT, ds_mem, ds_v);
    }
    catch (H5::Exception const& e) {
	throw exception("failed to read from HDF5 trajectory input file");
    }
}

} // namespace mdsim

#endif /* ! MDSIM_TRAJECTORY_HPP */
