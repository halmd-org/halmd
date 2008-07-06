/* Performance data
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

#ifndef MDSIM_PERF_HPP
#define MDSIM_PERF_HPP

#include <H5Cpp.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <sys/times.h>
#include <unistd.h>
#include "H5param.hpp"
#include "H5xx.hpp"
#include "accumulator.hpp"
#include "exception.hpp"
#include "log.hpp"
#include "statistics.hpp"


namespace mdsim
{

/**
 * performance class accumulators
 */
typedef boost::array<accumulator<double>, 5> perf_counters;

/**
 * performance data
 */
template <unsigned dimension, typename T>
class perf
{
public:
    perf() : m_offset(0), m_dirty(false) {}
    /** create HDF5 performance data output file */
    void open(std::string const& filename);
    /** returns HDF5 parameter group */
    H5param attrs();
    /** sample performance data */
    void sample(perf_counters const& times);
    /** clear performance counters */
    void commit();
    /** write performance data to HDF5 file */
    void flush(bool force);
    /** close HDF5 file */
    void close();

private:
    /** CPU tick accumulators */
    perf_counters m_times;
    /** HDF5 performance data output file */
    H5::H5File m_file;
    /** HDF5 datasets */
    boost::array<H5::DataSet, 5> m_dataset;
    /** HDF5 floating-point data type */
    H5::DataType m_tid;
    /** dataset offset */
    uint64_t m_offset;
    /** pending data bit */
    bool m_dirty;
};

/**
 * create HDF5 performance data output file
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::open(std::string const& filename)
{
    LOG("write performance data to file: " << filename);
    try {
	// truncate existing file
	m_file = H5::H5File(filename, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create performance data file");
    }
    // create parameter group
    m_file.createGroup("param");

    // floating-point data type
    m_tid = H5::PredType::NATIVE_DOUBLE;

    // extensible dataspace for performance data
    hsize_t dim[2] = { 0, 3 };
    hsize_t max_dim[2] = { H5S_UNLIMITED, 3 };
    hsize_t chunk_dim[2] = { 1, 3 };
    H5::DataSpace ds(2, dim, max_dim);
    H5::DSetCreatPropList cparms;
    cparms.setChunk(2, chunk_dim);

    try {
	H5::Group node(m_file.createGroup("times"));
	// CPU ticks for update cell lists
	m_dataset[0] = node.createDataSet("update_cells", m_tid, ds, cparms);
	// CPU ticks for update Verlet neighbour lists
	m_dataset[1] = node.createDataSet("update_neighbours", m_tid, ds, cparms);
	// CPU ticks for Lennard-Jones force update
	m_dataset[2] = node.createDataSet("update_forces", m_tid, ds, cparms);
	// CPU ticks for velocity-Verlet
	m_dataset[3] = node.createDataSet("velocity_verlet", m_tid, ds, cparms);
	// CPU ticks for MD simulation step
	m_dataset[4] = node.createDataSet("mdstep", m_tid, ds, cparms);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 performance datasets");
    }
}

/**
 * returns HDF5 parameter group
 */
template <unsigned dimension, typename T>
H5param perf<dimension, T>::attrs()
{
    return H5param(m_file.openGroup("param"));
}

/**
 * sample performance data
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::sample(perf_counters const& times)
{
    for (unsigned int i = 0; i < m_times.size(); ++i) {
	// accumulate values of accumulator
	m_times[i] += times[i];
    }
    m_dirty = true;
}

/**
 * output formatted accumulator times to stream
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, accumulator<T>& acc)
{
    os << std::setprecision(4) << (acc.mean() * 1000) << " ms";
    if (acc.count() > 1) {
	os << " (" << std::setprecision(4) << (acc.std() * 1000) << " ms, " << acc.count() << " calls)";
    }
    return os;
}

/**
 * commit HDF5 performance datasets
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::commit()
{
    LOG("mean CPU time for MD simulation step: " << m_times[4]);
    LOG("mean CPU time for velocity-Verlet integration: " << m_times[3]);
    LOG("mean CPU time for cell lists update: " << m_times[0]);
    LOG("mean CPU time for neighbour lists update: " << m_times[1]);
    LOG("mean CPU time for Lennard-Jones forces: " << m_times[2]);

    // write pending performance data to HDF5 file
    flush(false);

    for (unsigned int i = 0; i < m_times.size(); ++i) {
	// accumulate values of accumulator
	m_times[i].clear();
    }
    // increment offset in HDF5 datasets
    m_offset++;
    // clear pending data bit
    m_dirty = false;
}

/**
 * write performance data to HDF5 file
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::flush(bool force = true)
{
    if (!m_dirty)
	return;

    // file dataspace
    hsize_t dim[2] = { m_offset + 1, 3 };
    hsize_t count[2] = { 1, 1 };
    hsize_t start[2] = { m_offset, 0 };
    hsize_t stride[2] = { 1, 1 };
    hsize_t block[2] = { 1, 3 };
    H5::DataSpace ds(2, dim);
    ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    // memory dataspace
    H5::DataSpace mem_ds(2, block);
    boost::array<double, 3> data;

    try {
	for (unsigned int i = 0; i < m_times.size(); ++i) {
	    // write to HDF5 dataset
	    m_dataset[i].extend(dim);
	    data[0] = m_times[i].mean();
	    data[1] = m_times[i].std();
	    data[2] = m_times[i].count();
	    m_dataset[i].write(data.c_array(), m_tid, mem_ds, ds);
	}
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write performance data to HDF5 file");
    }

    if (force) {
	try {
	    m_file.flush(H5F_SCOPE_GLOBAL);
	}
	catch (H5::FileIException const& e) {
	    throw exception("failed to flush HDF5 performance file to disk");
	}
    }
}

/**
 * close HDF5 file
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::close()
{
    // write pending performance data to HDF5 file
    flush(false);

    try {
	m_file.close();
    }
    catch (H5::Exception const& e) {
	throw exception("failed to close performance data file");
    }
}

} // namespace mdsim

#endif /* ! MDSIM_PERFORMANCE_HPP */
