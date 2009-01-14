/* Performance data
 *
 * Copyright © 2008-2009  Peter Colberg
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

#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ljgpu/math/stat.hpp>
#include <ljgpu/sample/perf.hpp>
#include <ljgpu/util/exception.hpp>
#include <ljgpu/util/log.hpp>

namespace ljgpu
{

/**
 * performance counter descriptions
 */
perf::desc_map perf::desc = boost::assign::map_list_of
    ("boltzmann",		"Maxwell-Boltzmann distribution")
    ("event_queue",		"event queue processing")
    ("hilbert_sort",		"Hilbert space-filling curve sort")
    ("init_cells",		"cell lists initialisation")
    ("lattice",			"lattice generation")
    ("maximum_velocity",	"maximum velocity calculation")
    ("mdstep",			"MD simulation step")
    ("memcpy_cells",		"cell lists memcpy")
    ("potential_energy",	"potential energy sum calculation")
    ("reduce_squared_velocity",	"mean squared velocity calculation")
    ("reduce_velocity",		"velocity center of mass calculation")
    ("sample",			"phase space sampling")
    ("sample_memcpy",		"sample memcpy")
    ("update_cells",		"cell lists update")
    ("update_forces",		"Lennard-Jones force update")
    ("update_neighbours",	"neighbour lists update")
    ("velocity_verlet",		"velocity-Verlet integration")
    ("virial_sum",		"virial equation sum calculation")
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
    m_tid = H5::PredType::NATIVE_FLOAT;
    m_file.createGroup("param");
}

/**
 * sample performance data
 */
void perf::sample(counters const& times)
{
    BOOST_FOREACH(counter const& i, times) {
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
	os << " (" << std::fixed << std::setprecision(4)
	   << (acc.std() * 1000) << " ms, " << acc.count() << " calls)";
    }
    return os;
}

/**
 * commit HDF5 performance datasets
 */
void perf::commit()
{
    BOOST_FOREACH(counter const& i, m_times) {
	LOG(desc.at(i.first) << " average time: " << i.second);
    }

    // write pending performance data to HDF5 file
    flush(false);

    BOOST_FOREACH(counter& i, m_times) {
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

    if (m_offset == 0) {
	create_datasets();
    }

    write_datasets();

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

void perf::create_datasets()
{
    // extensible dataspace for performance data
    hsize_t dim[2] = { 0, 3 };
    hsize_t max_dim[2] = { H5S_UNLIMITED, 3 };
    hsize_t chunk_dim[2] = { 1, 3 };
    H5::DataSpace ds(2, dim, max_dim);
    H5::DSetCreatPropList cparms;
    cparms.setChunk(2, chunk_dim);

    try {
	H5::Group node(m_file.createGroup("times"));
	BOOST_FOREACH(counter const& i, m_times) {
	    std::string const& key = i.first;
	    if (m_dataset.find(key) == m_dataset.end()) {
		m_dataset[key] = node.createDataSet(key.c_str(), m_tid, ds, cparms);
	    }
	}
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 performance datasets");
    }
}

void perf::write_datasets()
{
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
	BOOST_FOREACH(counter const& i, m_times) {
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
}

} // namespace ljgpu
