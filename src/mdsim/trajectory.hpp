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
#include <boost/array.hpp>
#include <string>
#include "H5param.hpp"
#include "config.hpp"
#include "sample.hpp"

namespace mdsim {

template <bool writer>
class trajectory;

/**
 * trajectory file writer
 */
template <>
class trajectory<true>
{
public:
    /** create HDF5 trajectory output file */
    void open(std::string const& filename, unsigned int const& npart);
    /** close HDF5 trajectory output file */
    void close();
    /** flush HDF5 output file to disk */
    void flush();
    /** returns HDF5 parameter group */
    H5param attrs();
    /** write phase space sample */
    void sample(trajectory_sample const& sample, float_type const& time);

private:
    /** HDF5 trajectory output file */
    H5::H5File m_file;
    /** trajectory datasets for particle coordinates and velocities */
#ifdef USE_CUDA
    boost::array<H5::DataSet, 4> m_dataset;
#else
    boost::array<H5::DataSet, 3> m_dataset;
#endif
    /** memory dataspace for a single coordinates or velocities sample */
    H5::DataSpace m_ds_mem;
    /** file dataspace for a single coordinates or velocities sample */
    H5::DataSpace m_ds_file;
    /** file dataspace for simulation time */
    H5::DataSpace m_ds_scalar;
};

/**
 * trajectory file reader
 */
template <>
class trajectory<false>
{
public:
    /** open HDF5 trajectory input file */
    void open(std::string const& filename);
    /** close HDF5 trajectory input file */
    void close();
    /** read phase space sample */
    void read(std::vector<hvector>& r, std::vector<hvector>& v, int64_t index);

private:
    /** HDF5 trajectory input file */
    H5::H5File file;
};

} // namespace mdsim

#endif /* ! MDSIM_TRAJECTORY_HPP */
