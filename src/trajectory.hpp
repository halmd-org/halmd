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
#include "exception.hpp"
#include "options.hpp"


namespace mdsim {

template <typename T>
struct phase_space_point
{
    typedef T vector_type;

    /** coordinates of all particles in system */
    vector_type r;
    /** velocities of all particles in system */
    vector_type v;

    phase_space_point(uint64_t N) : r(N), v(N) { }
    phase_space_point() { }
};


template <unsigned dimension, typename T>
class trajectory
{
public:
    trajectory(options const& opts);
    void sample(phase_space_point<T> const& p, double const&, double const&);

    static void read(options const& opts, phase_space_point<T> &p);

private:
    H5::H5File file_;
    const uint64_t npart_;
    const uint64_t max_samples_;
    uint64_t samples_;
    boost::array<H5::DataSet, 2> dset_;
    H5::DataSpace ds_mem_;
    H5::DataSpace ds_file_;
};


/**
 * initialize HDF5 trajectory output file
 */
template <unsigned dimension, typename T>
trajectory<dimension, T>::trajectory(options const& opts) : npart_(opts.particles().value()), max_samples_(std::min(opts.steps().value(), opts.max_samples().value())), samples_(0)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    try {
	file_ = H5::H5File(opts.output_file_prefix().value() + ".trj", H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 trajectory file");
    }

    H5::DataSpace ds(H5S_SCALAR);
    H5::Group root(file_.openGroup("/"));

    unsigned int ndim = dimension;
    root.createAttribute("dimension", H5::PredType::NATIVE_UINT, ds).write(H5::PredType::NATIVE_UINT, &ndim);
    root.createAttribute("particles", H5::PredType::NATIVE_UINT64, ds).write(H5::PredType::NATIVE_UINT64, &npart_);
    root.createAttribute("steps", H5::PredType::NATIVE_UINT64, ds).write(H5::PredType::NATIVE_UINT64, &max_samples_);
    root.createAttribute("timestep", H5::PredType::NATIVE_DOUBLE, ds).write(H5::PredType::NATIVE_DOUBLE, &opts.timestep().value());
    // FIXME derived parameter box length is already calculated in ljfluid
    double box = pow(opts.particles().value() / opts.density().value(), 1.0 / dimension);
    root.createAttribute("box", H5::PredType::NATIVE_DOUBLE, ds).write(H5::PredType::NATIVE_DOUBLE, &box);

    hsize_t dim[3] = { max_samples_, npart_, dimension };
    ds_file_ = H5::DataSpace(3, dim);
    dset_[0] = file_.createDataSet("trajectory", H5::PredType::NATIVE_DOUBLE, ds_file_);
    dset_[1] = file_.createDataSet("velocity", H5::PredType::NATIVE_DOUBLE, ds_file_);

    hsize_t dim_mem[2] = { npart_, dimension };
    ds_mem_ = H5::DataSpace(2, dim_mem);
}

/**
 * write phase space sample to HDF5 dataset
 */
template <unsigned dimension, typename T>
void trajectory<dimension, T>::sample(phase_space_point<T> const& p, double const&, double const&)
{
    if (samples_ >= max_samples_)
	return;

    assert(p.r.size() == npart_);
    assert(p.v.size() == npart_);

    hsize_t count[3]  = { 1, npart_, 1 };
    hsize_t start[3]  = { samples_, 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, 1, dimension };

    ds_file_.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    // coordinates
    dset_[0].write(p.r.data(), H5::PredType::NATIVE_DOUBLE, ds_mem_, ds_file_);
    // velocities
    dset_[1].write(p.v.data(), H5::PredType::NATIVE_DOUBLE, ds_mem_, ds_file_);

    samples_++;
}

/**
 * read phase space sample from HDF5 trajectory input file
 */
template <unsigned dimension, typename T>
void trajectory<dimension, T>::read(options const& opts, phase_space_point<T> &p)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    H5::H5File file;

    try {
	file = H5::H5File(opts.trajectory_input_file().value(), H5F_ACC_RDONLY);
    }
    catch (H5::Exception const& e) {
	throw exception("failed to open HDF5 trajectory input file");
    }

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
	hsize_t nsamples = dim_r[0];
	// number of particles
	hsize_t npart = dim_r[1];
	// positional coordinate dimension
	hsize_t ndim = dim_r[2];

	// validate dataspace dimensions
	if (nsamples < 1) {
	    throw exception("trajectory input file has invalid number of samples");
	}
	if (npart < 1) {
	    throw exception("trajectory input file has invalid number of particles");
	}
	if (ndim != dimension) {
	    throw exception("trajectory input file has invalid coordinate dimension");
	}

	// check if sample number is within bounds
	if ((opts.sample().value() >= int(nsamples)) || ((-opts.sample().value()) > int(nsamples))) {
	    throw exception("trajectory input sample number out of bounds");
	}
	hsize_t sample = (opts.sample().value() < 0) ? (opts.sample().value() + int(nsamples)) : opts.sample().value();

	// allocate memory for sample
	p.r.resize(npart);
	p.v.resize(npart);

	// read sample from dataset
	hsize_t dim_mem[2] = { npart, ndim };
	H5::DataSpace ds_mem(2, dim_mem);

	hsize_t count[3]  = { 1, npart, 1 };
	hsize_t start[3]  = { sample, 0, 0 };
	hsize_t stride[3] = { 1, 1, 1 };
	hsize_t block[3]  = { 1, 1, ndim };

	// coordinates
	ds_r.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dset_r.read(p.r.data(), H5::PredType::NATIVE_DOUBLE, ds_mem, ds_r);
	// velocities
	ds_v.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	dset_v.read(p.v.data(), H5::PredType::NATIVE_DOUBLE, ds_mem, ds_v);
    }
    catch (H5::Exception const& e) {
	throw exception("failed to read from HDF5 trajectory input file");
    }
}

} // namespace mdsim

#endif /* ! MDSIM_TRAJECTORY_HPP */
