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
#include "accumulator.hpp"
#include "ljfluid.hpp"
#include "options.hpp"
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
    energy(options const& opts);
    /** write global simulation parameters to thermodynamic equilibrium properties output file */
    void write_param(H5param const& param);
    void sample(phase_space_point<std::vector<T> > const& p, double const& en_pot, double const& virial, uint64_t);
    void write();

private:
    /** simulation timestep */
    double timestep_;
    /** particle density */
    double density_;
    /** sample count */
    unsigned int samples_;
    /** number of samples */
    unsigned int max_samples_;

    /** thermodynamic equilibrium properties output file */
    H5::H5File file_;

    /** thermodynamic equilibrium properties */
    std::vector<double> en_pot_;
    std::vector<double> en_kin_;
    std::vector<double> en_tot_;
    std::vector<double> temp_;
    std::vector<double> press_;
    std::vector<T> v_cm_;
};


template <unsigned dimension, typename T>
energy<dimension, T>::energy(options const& opts) : timestep_(opts.timestep().value()), density_(opts.density().value()), samples_(0)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    try {
	file_ = H5::H5File(opts.output_file_prefix().value() + ".tep", H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 energy file");
    }

    // number of samples
    max_samples_ = std::min(opts.max_samples().value(), opts.steps().value());

    try {
	en_pot_.reserve(max_samples_);
	en_kin_.reserve(max_samples_);
	en_tot_.reserve(max_samples_);
	temp_.reserve(max_samples_);
	press_.reserve(max_samples_);
	v_cm_.reserve(max_samples_);
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate thermodynamic equilibrium properties buffer");
    }
}

/**
 * write global simulation parameters to thermodynamic equilibrium properties output file
 */
template <unsigned dimension, typename T>
void energy<dimension, T>::write_param(H5param const& param)
{
    param.write(file_.createGroup("/parameters"));
}
 
/**
 * sample thermodynamic equilibrium properties
 */
template <unsigned dimension, typename T>
void energy<dimension, T>::sample(phase_space_point<std::vector<T> > const& p, double const& en_pot, double const& virial, uint64_t)
{
    if (samples_ >= max_samples_) return;

    // mean squared velocity
    accumulator<double> vv;
    foreach (T const& v, p.v) {
	vv += v * v;
    }

    // mean potential energy per particle
    en_pot_.push_back(en_pot);
    // mean kinetic energy per particle
    en_kin_.push_back(vv.mean() / 2);
    // mean total energy per particle
    en_tot_.push_back(en_pot_.back() + en_kin_.back());
    // temperature
    temp_.push_back(vv.mean() / dimension);
    // pressure
    press_.push_back(density_ * (vv.mean() + virial));
    // velocity center of mass
    v_cm_.push_back(mean(p.v.begin(), p.v.end()));

    samples_++;
}


/**
 * write thermodynamic equilibrium properties buffer to file
 */
template <unsigned dimension, typename T>
void energy<dimension, T>::write()
{
    // create dataspaces for scalar and vector types
    hsize_t dim_scalar[2] = { max_samples_, 1 };
    hsize_t dim_vector[2] = { max_samples_, dimension };
    H5::DataSpace ds_scalar(2, dim_scalar);
    H5::DataSpace ds_vector(2, dim_vector);

    // HDF5 datasets for thermodynamic equilibrium properties
    boost::array<H5::DataSet, 6> dataset_;
    H5::DataType tid(H5::PredType::NATIVE_DOUBLE);

    try {
	// mean potential energy per particle
	dataset_[0] = file_.createDataSet("EPOT", tid, ds_scalar);
	// mean kinetic energy per particle
	dataset_[1] = file_.createDataSet("EKIN", tid, ds_scalar);
	// mean total energy per particle
	dataset_[2] = file_.createDataSet("ETOT", tid, ds_scalar);
	// temperature
	dataset_[3] = file_.createDataSet("TEMP", tid, ds_scalar);
	// pressure
	dataset_[4] = file_.createDataSet("PRESS", tid, ds_scalar);
	// velocity center of mass
	dataset_[5] = file_.createDataSet("VCM", tid, ds_vector);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create datasets in HDF5 energy file");
    }

    try {
	// mean potential energy per particle
	dataset_[0].write(en_pot_.data(), tid, ds_scalar, ds_scalar);
	// mean kinetic energy per particle
	dataset_[1].write(en_kin_.data(), tid, ds_scalar, ds_scalar);
	// mean total energy per particle
	dataset_[2].write(en_tot_.data(), tid, ds_scalar, ds_scalar);
	// temperature
	dataset_[3].write(temp_.data(), tid, ds_scalar, ds_scalar);
	// pressure
	dataset_[4].write(press_.data(), tid, ds_scalar, ds_scalar);
	// velocity center of mass
	dataset_[5].write(v_cm_.data(), tid, ds_vector, ds_vector);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write thermodynamic equilibrium properties to HDF5 energy file");
    }
}

} // namespace mdsim

#undef foreach

#endif /* ! MDSIM_ENERGY_HPP */
