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

public:
    energy() : m_samples(0) {}
    /** create HDF5 thermodynamic equilibrium properties output file */
    void open(std::string const& filename);
    /** returns HDF5 parameter group */
    H5param attrs();
    /** sample thermodynamic equilibrium properties */
    void sample(std::vector<T> const& v, double const& virial, double const& density, double const& timestep, double const& time);
    /** write thermodynamic equilibrium properties to HDF5 file */
    void write();
    /** close HDF5 thermodynamic equilibrium properties output file */
    void close();

private:
    /** number of aquired samples */
    uint64_t m_samples;

    /** thermodynamic equilibrium properties */
    std::vector<scalar_pair> m_en_kin;
    std::vector<scalar_pair> m_temp;
    std::vector<scalar_pair> m_press;
    std::vector<vector_pair> m_v_cm;
    /** HDF5 thermodynamic equilibrium properties output file */
    H5::H5File m_file;
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
}


/**
 * write thermodynamic equilibrium properties to HDF5 file
 */
template <unsigned dimension, typename T>
void energy<dimension, T>::write()
{
    // create dataspaces for scalar and vector types
    hsize_t dim_scalar[2] = { m_samples, 2 };
    hsize_t dim_vector[2] = { m_samples, 1 + dimension };
    H5::DataSpace ds_scalar(2, dim_scalar);
    H5::DataSpace ds_vector(2, dim_vector);

    // HDF5 datasets for thermodynamic equilibrium properties
    boost::array<H5::DataSet, 4> dataset_;
    H5::DataType tid(H5::PredType::NATIVE_DOUBLE);

    try {
	// mean kinetic energy per particle
	dataset_[0] = m_file.createDataSet("EKIN", tid, ds_scalar);
	// temperature
	dataset_[1] = m_file.createDataSet("TEMP", tid, ds_scalar);
	// pressure
	dataset_[2] = m_file.createDataSet("PRESS", tid, ds_scalar);
	// velocity center of mass
	dataset_[3] = m_file.createDataSet("VCM", tid, ds_vector);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create datasets in HDF5 energy file");
    }

    try {
	// mean kinetic energy per particle
	dataset_[0].write(m_en_kin.data(), tid, ds_scalar, ds_scalar);
	// temperature
	dataset_[1].write(m_temp.data(), tid, ds_scalar, ds_scalar);
	// pressure
	dataset_[2].write(m_press.data(), tid, ds_scalar, ds_scalar);
	// velocity center of mass
	dataset_[3].write(m_v_cm.data(), tid, ds_vector, ds_vector);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write thermodynamic equilibrium properties to HDF5 energy file");
    }
}

/**
 * close HDF5 file
 */
template <unsigned dimension, typename T>
void energy<dimension, T>::close()
{
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
