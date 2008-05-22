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
};


template <unsigned dimension, typename T>
class trajectory
{
public:
    trajectory(options const& opts);
    void sample(phase_space_point<T> const& p, double const&, double const&);

private:
    H5::H5File file_;
    const uint64_t npart_;
    const uint64_t max_samples_;
    uint64_t samples_;
    H5::DataSpace ds_;
    H5::DataSet dset_[2];
    H5::DataSpace ds_src_;
    H5::DataSpace ds_dst_;
};


template <unsigned dimension, typename T>
trajectory<dimension, T>::trajectory(options const& opts) : npart_(opts.npart()), max_samples_(std::min(opts.steps(), opts.max_samples())), samples_(0)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    try {
	file_ = H5::H5File(opts.trajectory_output_file(), H5F_ACC_TRUNC);
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
    root.createAttribute("timestep", H5::PredType::NATIVE_DOUBLE, ds).write(H5::PredType::NATIVE_DOUBLE, &opts.timestep());
    // FIXME derived parameter box length is already calculated in ljfluid
    double box = pow(opts.npart() / opts.density(), 1.0 / dimension);
    root.createAttribute("box", H5::PredType::NATIVE_DOUBLE, ds).write(H5::PredType::NATIVE_DOUBLE, &box);

    hsize_t dim1[3] = { max_samples_, npart_, dimension };
    ds_ = H5::DataSpace(3, dim1);
    dset_[0] = file_.createDataSet("trajectory", H5::PredType::NATIVE_DOUBLE, ds_);
    dset_[1] = file_.createDataSet("velocity", H5::PredType::NATIVE_DOUBLE, ds_);

    hsize_t dim2[2] = { npart_, dimension };
    ds_src_ = H5::DataSpace(2, dim2);
    ds_dst_ = ds_;
}

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

    ds_dst_.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    // coordinates
    dset_[0].write(p.r.data(), H5::PredType::NATIVE_DOUBLE, ds_src_, ds_dst_);
    // velocities
    dset_[1].write(p.v.data(), H5::PredType::NATIVE_DOUBLE, ds_src_, ds_dst_);

    samples_++;
}

} // namespace mdsim

#endif /* ! MDSIM_TRAJECTORY_HPP */
