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

#include "exception.hpp"
#include "log.hpp"
#include "perf.hpp"
#include "statistics.hpp"
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <iomanip>
#include <iostream>
#include <limits>

namespace mdsim
{

/**
 * performance counter datasets
 */
typedef boost::unordered_map<std::string, std::string> perf_dataset_map;
static perf_dataset_map const perf_dataset = boost::assign::map_list_of
#if USE_CUDA
    ("mdstep",			"GPU time for MD simulation step")
    ("velocity_verlet",		"GPU time for velocity-Verlet integration")
# if defined(USE_NEIGHBOUR) && defined(USE_HILBERT_ORDER)
    ("hilbert_sort",		"GPU time for Hilbert space-filling curve sort")
# endif
# if defined(USE_NEIGHBOUR) || defined(USE_CELL)
    ("update_cells",		"GPU time for cell lists update")
# endif
# if defined(USE_NEIGHBOUR)
    ("update_neighbours",	"GPU time for neighbour lists update")
    ("maximum_velocity",	"GPU time for maximum velocity calculation")
# endif
# if defined(USE_CELL)
    ("init_cells",		"GPU time for cell lists initialisation")
    ("memcpy_cells",		"GPU time for cell lists memcpy")
# endif
    ("update_forces",		"GPU time for Lennard-Jones force update")
# if defined(USE_NEIGHBOUR) || !defined(USE_CELL)
    ("potential_energy",	"GPU time for potential energy sum calculation")
    ("virial_sum",		"GPU time for virial equation sum calculation")
# endif
    ("sample_memcpy",		"GPU time for sample memcpy")
    ("lattice",			"GPU time for lattice generation")
    ("boltzmann",		"GPU time for Maxwell-Boltzmann distribution")
#else /* ! USE_CUDA */
    ("update_cells",		"CPU time for cell lists update")
    ("update_neighbours",	"CPU time for neighbour lists update")
    ("update_forces",		"CPU time for Lennard-Jones force update")
    ("velocity_verlet",		"CPU time for velocity-Verlet integration")
    ("mdstep",			"CPU time for MD simulation step")
#endif /* ! USE_CUDA */
    ;

/**
 * create HDF5 performance data output file
 */
void perf::open(std::string const& filename)
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
	// FIXME only create necessary (implementation-dependent) datasets
	BOOST_FOREACH(perf_dataset_map::value_type const& i, perf_dataset) {
	    m_dataset[i.first] = node.createDataSet(i.first.c_str(), m_tid, ds, cparms);
	}
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 performance datasets");
    }
}

/**
 * returns HDF5 parameter group
 */
H5param perf::attrs()
{
    return H5param(m_file.openGroup("param"));
}

/**
 * sample performance data
 */
void perf::sample(perf_counters const& times)
{
    BOOST_FOREACH(perf_counters::value_type const& i, times) {
	// accumulate values of accumulator
	m_times[i.first] += i.second;
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
void perf::commit()
{
    BOOST_FOREACH(perf_counters::value_type const& i, m_times) {
	LOG("mean " << perf_dataset.at(i.first) << ": " << i.second);
    }

    // write pending performance data to HDF5 file
    flush(false);

    BOOST_FOREACH(perf_counters::value_type& i, m_times) {
	// accumulate values of accumulator
	i.second.clear();
    }
    // increment offset in HDF5 datasets
    m_offset++;
    // clear pending data bit
    m_dirty = false;
}

/**
 * write performance data to HDF5 file
 */
void perf::flush(bool force)
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
	BOOST_FOREACH(perf_counters::value_type const& i, m_times) {
	    // write to HDF5 dataset
	    m_dataset[i.first].extend(dim);
	    data[0] = i.second.mean();
	    data[1] = i.second.std();
	    data[2] = i.second.count();
	    m_dataset[i.first].write(data.c_array(), m_tid, mem_ds, ds);
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
void perf::close()
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
