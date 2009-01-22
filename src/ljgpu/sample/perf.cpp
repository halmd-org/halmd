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

#include <algorithm>
#include <boost/array.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ljgpu/math/stat.hpp>
#include <ljgpu/sample/perf.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <ljgpu/util/exception.hpp>
#include <ljgpu/util/log.hpp>

#define foreach BOOST_FOREACH

namespace ljgpu
{

/**
 * performance counter descriptions
 */
perf::desc_map perf::m_desc = boost::assign::map_list_of
    ("anderson_thermostat",	"Anderson thermostat")
    ("boltzmann",		"Boltzmann distribution")
    ("event_queue",		"event queue processing")
    ("hilbert_sort",		"Hilbert curve sort")
    ("init_cells",		"cell lists initialisation")
    ("lattice",			"lattice generation")
    ("maximum_velocity",	"maximum velocity reduction")
    ("mdstep",			"MD simulation step")
    ("memcpy_cells",		"cell lists memcpy")
    ("potential_energy",	"potential energy sum reduction")
    ("reduce_squared_velocity",	"mean squared velocity reduction")
    ("reduce_velocity",		"velocity center of mass reduction")
    ("sample",			"phase space sampling")
    ("sample_memcpy",		"sample memcpy")
    ("update_cells",		"cell lists update")
    ("update_forces",		"Lennard-Jones force update")
    ("update_neighbours",	"neighbour lists update")
    ("velocity_verlet",		"velocity-Verlet integration")
    ("virial_sum",		"virial equation sum reduction")
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
}

/**
 * sample performance data
 */
void perf::sample(counters const& times)
{
    foreach (counter const& i, times) {
	// accumulate values of accumulator
	m_times[i.first] += i.second;
    }
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
 * output formatted performance statistics to stream
 */
std::ostream& operator<<(std::ostream& os, perf::counters const& times)
{
    std::vector<std::string> keys;
    foreach (perf::counter const& i, times) {
	keys.push_back(i.first);
    }
    std::sort(keys.begin(), keys.end());
    foreach (std::string const& key, keys) {
	os << perf::desc(key) << ": " << times.at(key) << std::endl;
    }
    return os;
}

/**
 * write performance data to HDF5 file
 */
void perf::flush()
{
    H5::Group node;
    try {
	H5XX_NO_AUTO_PRINT(H5::FileIException);
	node = m_file.openGroup("times");
    }
    catch (H5::FileIException const&) {
	node = m_file.createGroup("times");
    }

    // dataspace for performance data
    hsize_t dim[2] = { 1, 3 };
    H5::DataSpace ds_file(2, dim);

    try {
	foreach (counter const& i, m_times) {
	    std::string const& key = i.first;
	    if (m_dataset.find(key) == m_dataset.end()) {
		m_dataset[key] = node.createDataSet(key.c_str(), m_tid, ds_file);
	    }
	}
    }
    catch (H5::FileIException const&) {
	throw exception("failed to create HDF5 performance datasets");
    }

    // memory dataspace
    H5::DataSpace ds_mem(2, dim);
    boost::array<float, 3> data;
    try {
	foreach (counter const& i, m_times) {
	    // write to HDF5 dataset
	    m_dataset[i.first].extend(dim);
	    data[0] = i.second.mean();
	    data[1] = i.second.std();
	    data[2] = i.second.count();
	    m_dataset[i.first].write(data.c_array(), m_tid, ds_mem, ds_file);
	}
    }
    catch (H5::FileIException const&) {
	throw exception("failed to write performance data to HDF5 file");
    }

    try {
	m_file.flush(H5F_SCOPE_GLOBAL);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to flush HDF5 performance file to disk");
    }
}

/**
 * close HDF5 file
 */
void perf::close()
{
    // write pending performance data to HDF5 file
    flush();

    try {
	m_file.close();
    }
    catch (H5::Exception const& e) {
	throw exception("failed to close performance data file");
    }
}

} // namespace ljgpu
