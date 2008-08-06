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
#include "H5param.hpp"
#include "H5xx.hpp"
#include "accumulator.hpp"
#include "exception.hpp"
#include "log.hpp"
#include "statistics.hpp"

namespace mdsim
{

/**
 * performance counter types
 */
enum {
    GPU_TIME_MDSTEP,
    GPU_TIME_VELOCITY_VERLET,
#ifdef USE_CELL
    GPU_TIME_HILBERT_SORT,
    GPU_TIME_UPDATE_CELLS,
    GPU_TIME_UPDATE_NEIGHBOURS,
#endif
    GPU_TIME_UPDATE_FORCES,
    GPU_TIME_SAMPLE_MEMCPY,
    GPU_TIME_LATTICE,
    GPU_TIME_BOLTZMANN,
    NUM_PERF_COUNTERS,
};

/**
 * performance counter datasets
 */
static char const* perf_dataset[][2] = {
    { "mdstep",			"GPU time for MD simulation step" },
    { "velocity_verlet",	"GPU time for velocity-Verlet integration" },
#ifdef USE_CELL
    { "hilbert_sort",		"GPU time for Hilbert space-filling curve sort" },
    { "update_cells",		"GPU time for cell lists update" },
    { "update_neighbours",	"GPU time for neighbour lists update" },
#endif
    { "update_forces",		"GPU time for Lennard-Jones force update" },
    { "sample_memcpy",		"GPU time for sample memcpy" },
    { "lattice",		"GPU time for lattice generation" },
    { "boltzmann",		"GPU time for Maxwell-Boltzmann distribution" },
};

/**
 * performance class accumulators
 */
typedef boost::array<accumulator<float>, NUM_PERF_COUNTERS> perf_counters;

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
    boost::array<H5::DataSet, perf_counters::static_size> m_dataset;
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
    m_tid = H5::PredType::NATIVE_FLOAT;

    // extensible dataspace for performance data
    hsize_t dim[2] = { 0, 3 };
    hsize_t max_dim[2] = { H5S_UNLIMITED, 3 };
    hsize_t chunk_dim[2] = { 1, 3 };
    H5::DataSpace ds(2, dim, max_dim);
    H5::DSetCreatPropList cparms;
    cparms.setChunk(2, chunk_dim);

    try {
	H5::Group node(m_file.createGroup("times"));
	for (unsigned int i = 0; i < m_dataset.size(); ++i) {
	    m_dataset[i] = node.createDataSet(perf_dataset[i][0], m_tid, ds, cparms);
	}
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
 * output formatted accumulator values to stream
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, accumulator<T> const& acc)
{
    os << std::fixed << std::setprecision(4) << (acc.mean() * 1000) << " ms";
    if (acc.count() > 1) {
	os << " (" << std::fixed << std::setprecision(4) << (acc.std() * 1000) << " ms, " << acc.count() << " calls)";
    }
    return os;
}

/**
 * commit HDF5 performance datasets
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::commit()
{
    for (unsigned int i = 0; i < m_times.size(); ++i) {
	LOG("mean " << perf_dataset[i][1] << ": " << m_times[i]);
    }

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
    boost::array<float, 3> data;
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

#endif /* ! MDSIM_PERF_HPP */
