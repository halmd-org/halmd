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
#include "exception.hpp"
#include <assert.h>


namespace mdsim {

template <unsigned int NDIM, typename T>
class trajectory
{
public:
    trajectory(char const* path, unsigned int npart, unsigned int steps);

    void write(T const& coord, T const& vel);

private:
    H5::H5File file_;
    unsigned int sets_;
    const unsigned int npart_;
    const unsigned int steps_;
    H5::DataSpace ds_;
    H5::DataSet dset_[2];
    H5::DataSpace ds_src_;
    H5::DataSpace ds_dst_;
};


template <unsigned int NDIM, typename T>
trajectory<NDIM, T>::trajectory(char const* path, unsigned int npart, unsigned int steps) : sets_(0), npart_(npart), steps_(steps)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    try {
	file_ = H5::H5File(path, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 trajectory file");
    }

    hsize_t dim[1] = { 1 };
    H5::DataSpace ds(1, dim);
    H5::FloatType dt(H5::PredType::NATIVE_UINT);
    H5::Group root(file_.openGroup("/"));

    unsigned int ndim = NDIM;
    root.createAttribute("dimensions", dt, ds).write(dt, &ndim);
    root.createAttribute("particles", dt, ds).write(dt, &npart_);
    root.createAttribute("steps", dt, ds).write(dt, &steps_);

    hsize_t dim1[3] = { steps_, npart_, NDIM };
    ds_ = H5::DataSpace(3, dim1);
    dset_[0] = file_.createDataSet("trajectory", H5::PredType::NATIVE_DOUBLE, ds_);
    dset_[1] = file_.createDataSet("velocity", H5::PredType::NATIVE_DOUBLE, ds_);

    hsize_t dim2[2] = { npart_, NDIM };
    ds_src_ = H5::DataSpace(2, dim2);
    ds_dst_ = ds_;
}

template <unsigned int NDIM, typename T>
void trajectory<NDIM, T>::write(T const& coord, T const& vel)
{
    assert(sets_ < steps_);
    assert(coord.size() == npart_);
    assert(vel.size() == npart_);

    hsize_t count[3]  = { 1, npart_, 1 };
    hsize_t start[3]  = { sets_, 0, 0 };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, 1, NDIM };

    ds_dst_.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    // coordinates
    dset_[0].write(coord.data(), H5::PredType::NATIVE_DOUBLE, ds_src_, ds_dst_);
    // velocities
    dset_[1].write(vel.data(), H5::PredType::NATIVE_DOUBLE, ds_src_, ds_dst_);

    sets_++;
}

} // namespace mdsim

#endif /* ! MDSIM_TRAJECTORY_HPP */
