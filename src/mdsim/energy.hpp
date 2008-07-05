/* Thermodynamic equilibrium properties
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

#ifndef MDSIM_ENERGY_HPP
#define MDSIM_ENERGY_HPP

#include <H5Cpp.h>
#include <algorithm>
#include <boost/foreach.hpp>
#include <string>
#include <vector>
#include "H5param.hpp"
#include "accumulator.hpp"
#include "log.hpp"
#include "statistics.hpp"

#define foreach BOOST_FOREACH

namespace mdsim
{

/**
 * Thermodynamic equilibrium properties
 */
template <unsigned dimension, typename T>
class energy
{
public:
    /** time and scalar property */
    typedef std::pair<double, double> scalar_pair;
    /** time and vector property */
    typedef std::pair<double, T> vector_pair;

    enum {
	/** HDF5 dataset chunk size */
	CHUNK_SIZE = 2500,
    };

public:
    energy() : m_samples(0), m_samples_buffer(0), m_samples_file(0) {}
    /** create HDF5 thermodynamic equilibrium properties output file */
    void open(std::string const& filename);
    /** returns HDF5 parameter group */
    H5param attrs();
    /** sample thermodynamic equilibrium properties */
    void sample(std::vector<T> const& v, double const& virial, double const& density, double const& timestep, double const& time);
    /** close HDF5 thermodynamic equilibrium properties output file */
    void close();

private:
    /** write thermodynamic equilibrium properties to HDF5 file */
    void flush();

private:
    /** number of aquired samples */
    uint64_t m_samples;
    /** number of samples in memory */
    uint64_t m_samples_buffer;
    /** number of samples committed to file */
    uint64_t m_samples_file;

    /** thermodynamic equilibrium properties */
    std::vector<scalar_pair> m_en_kin;
    std::vector<scalar_pair> m_temp;
    std::vector<scalar_pair> m_press;
    std::vector<vector_pair> m_v_cm;
    /** HDF5 thermodynamic equilibrium properties output file */
    H5::H5File m_file;
    /** HDF5 datasets */
    boost::array<H5::DataSet, 4> m_dataset;
    /** HDF5 floating-point data type */
    H5::DataType m_tid;
};

/**
 * create HDF5 thermodynamic equilibrium properties output file
 */
template <unsigned dimension, typename T>
void energy<dimension, T>::open(std::string const& filename)
{
    LOG("write thermodynamic equilibrium properties to file: " << filename);
    try {
	// truncate existing file
	m_file = H5::H5File(filename, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create thermodynamic equilibrium properties output file");
    }
    // create parameter group
    m_file.createGroup("param");

    // extensible dataspace for scalar properties
    hsize_t scalar_dim[2] = { 0, 2 };
    hsize_t scalar_max_dim[2] = { H5S_UNLIMITED, 2 };
    hsize_t scalar_chunk_dim[2] = { CHUNK_SIZE, 2 };
    H5::DataSpace scalar_ds(2, scalar_dim, scalar_max_dim);
    H5::DSetCreatPropList scalar_cparms;
    scalar_cparms.setChunk(2, scalar_chunk_dim);

    // extensible dataspace for vector properties
    hsize_t vector_dim[2] = { 0, dimension + 1 };
    hsize_t vector_max_dim[2] = { H5S_UNLIMITED, dimension + 1 };
    hsize_t vector_chunk_dim[2] = { CHUNK_SIZE, dimension + 1 };
    H5::DataSpace vector_ds(2, vector_dim, vector_max_dim);
    H5::DSetCreatPropList vector_cparms;
    vector_cparms.setChunk(2, vector_chunk_dim);

    // floating-point data type
    m_tid = H5::PredType::NATIVE_DOUBLE;

    try {
	// mean kinetic energy per particle
	m_dataset[0] = m_file.createDataSet("EKIN", m_tid, scalar_ds, scalar_cparms);
	// temperature
	m_dataset[1] = m_file.createDataSet("TEMP", m_tid, scalar_ds, scalar_cparms);
	// pressure
	m_dataset[2] = m_file.createDataSet("PRESS", m_tid, scalar_ds, scalar_cparms);
	// velocity center of mass
	m_dataset[3] = m_file.createDataSet("VCM", m_tid, vector_ds, vector_cparms);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create datasets in HDF5 energy file");
    }

    try {
	// allocate thermodynamic equilibrium property buffers
	m_en_kin.reserve(CHUNK_SIZE);
	m_temp.reserve(CHUNK_SIZE);
	m_press.reserve(CHUNK_SIZE);
	m_v_cm.reserve(CHUNK_SIZE);
    }
    catch (std::bad_alloc const&) {
	throw exception("failed to allocate thermodynamic equilibrium property buffers");
    }
}

/**
 * returns HDF5 parameter group
 */
template <unsigned dimension, typename T>
H5param energy<dimension, T>::attrs()
{
    return H5param(m_file.openGroup("param"));
}

/**
 * sample thermodynamic equilibrium properties
 */
template <unsigned dimension, typename T>
void energy<dimension, T>::sample(std::vector<T> const& v, double const& virial, double const& density, double const& timestep, double const& time)
{
    // mean squared velocity
    accumulator<double> vv;
    foreach (T const& v_, v) {
	vv += v_ * v_;
    }

    // mean kinetic energy per particle
    m_en_kin.push_back(scalar_pair(time, vv.mean() / 2));
    // temperature
    m_temp.push_back(scalar_pair(time, vv.mean() / dimension));
    // pressure
    m_press.push_back(scalar_pair(time, density / dimension * (vv.mean() + virial / timestep)));
    // velocity center of mass
    m_v_cm.push_back(vector_pair(time, mean(v.begin(), v.end())));

    m_samples++;
    m_samples_buffer++;

    // commit full buffers to file
    if (m_samples_buffer >= CHUNK_SIZE) {
	flush();
    }
}

/**
 * write thermodynamic equilibrium properties to HDF5 file
 */
template <unsigned dimension, typename T>
void energy<dimension, T>::flush()
{
    // file dataspaces
    hsize_t scalar_dim[2] = { m_samples, 2 };
    hsize_t vector_dim[2] = { m_samples, dimension + 1 };
    H5::DataSpace scalar_ds(2, scalar_dim);
    H5::DataSpace vector_ds(2, vector_dim);

    // file dataspace hyperslabs
    hsize_t count[2] = { 1, 1 };
    hsize_t start[2] = { m_samples_file, 0 };
    hsize_t stride[2] = { 1, 1 };
    hsize_t scalar_block[2] = { m_samples_buffer, 2 };
    hsize_t vector_block[2] = { m_samples_buffer, dimension + 1 };
    scalar_ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, scalar_block);
    vector_ds.selectHyperslab(H5S_SELECT_SET, count, start, stride, vector_block);

    // memory dataspaces
    H5::DataSpace scalar_mem_ds(2, scalar_block);
    H5::DataSpace vector_mem_ds(2, vector_block);

    try {
	// mean kinetic energy per particle
	m_dataset[0].extend(scalar_dim);
	m_dataset[0].write(m_en_kin.data(), m_tid, scalar_mem_ds, scalar_ds);
	m_en_kin.clear();
	// temperature
	m_dataset[1].extend(scalar_dim);
	m_dataset[1].write(m_temp.data(), m_tid, scalar_mem_ds, scalar_ds);
	m_temp.clear();
	// pressure
	m_dataset[2].extend(scalar_dim);
	m_dataset[2].write(m_press.data(), m_tid, scalar_mem_ds, scalar_ds);
	m_press.clear();
	// velocity center of mass
	m_dataset[3].extend(vector_dim);
	m_dataset[3].write(m_v_cm.data(), m_tid, vector_mem_ds, vector_ds);
	m_v_cm.clear();
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write thermodynamic equilibrium properties to HDF5 energy file");
    }

    m_samples_file = m_samples;
    m_samples_buffer = 0;
}

/**
 * close HDF5 file
 */
template <unsigned dimension, typename T>
void energy<dimension, T>::close()
{
    // commit remaining samples to file
    if (m_samples_buffer > 0) {
	flush();
    }

    try {
	m_file.close();
    }
    catch (H5::Exception const& e) {
	throw exception("failed to close HDF5 correlations output file");
    }
}

} // namespace mdsim

#undef foreach

#endif /* ! MDSIM_ENERGY_HPP */
