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

#include <stdint.h>
#include <vector>
#include "accumulator.hpp"
#include "ljfluid.hpp"
#include "options.hpp"
#include "statistics.hpp"


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
    void sample(phase_space_point<T> const& p, double const& en_pot, double const& virial);

private:
    /** HDF5 energy file */
    H5::H5File file_;
    /** simulation timestep */
    double timestep_;
    /** particle density */
    double density_;
    /** number of samples */
    uint64_t samples_;

    /** HDF5 data sets for thermodynamic equilibrium properties */
    std::vector<H5::DataSet> dset_;
    /** HDF5 data space for selection within data set */
    H5::DataSpace ds_scalar_, ds_vector_;
};


template <unsigned dimension, typename T>
energy<dimension, T>::energy(options const& opts) : timestep_(opts.timestep()), density_(opts.density()), samples_(0)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    try {
	file_ = H5::H5File(opts.energy_output_file(), H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 energy file");
    }

    H5::DataType dt(H5::PredType::NATIVE_DOUBLE);
    H5::Group root(file_.openGroup("/"));

    root.createAttribute("timestep", dt, H5S_SCALAR).write(dt, &timestep_);

    hsize_t dim_scalar[2] = { opts.steps(), 1 };
    hsize_t dim_vector[2] = { opts.steps(), dimension };
    ds_scalar_ = H5::DataSpace(2, dim_scalar);
    ds_vector_ = H5::DataSpace(2, dim_vector);

    // mean potential energy per particle
    dset_.push_back(file_.createDataSet("EPOT", dt, ds_scalar_));
    // mean kinetic energy per particle
    dset_.push_back(file_.createDataSet("EKIN", dt, ds_scalar_));
    // mean total energy per particle
    dset_.push_back(file_.createDataSet("ETOT", dt, ds_scalar_));
    // temperature
    dset_.push_back(file_.createDataSet("TEMP", dt, ds_scalar_));
    // pressure
    dset_.push_back(file_.createDataSet("PRESS", dt, ds_scalar_));
    // velocity center of mass
    dset_.push_back(file_.createDataSet("VCM", dt, ds_vector_));
}

template <unsigned dimension, typename T>
void energy<dimension, T>::sample(phase_space_point<T> const& p, double const& en_pot, double const& virial)
{
    // mean squared velocity
    accumulator<double> vv;
    for (typename T::const_iterator it = p.v.begin(); it != p.v.end(); ++it) {
	vv += *it * *it;
    }

    // mean kinetic energy per particle
    double en_kin = vv.mean() / 2.;
    // mean total energy per particle
    double en_tot = en_pot + en_kin;
    // temperature
    double temp = vv.mean() / dimension;
    // pressure
    double press = density_ * (vv.mean() + virial);
    // velocity center of mass
    typename T::value_type v_cm = mean(p.v.begin(), p.v.end());

    hsize_t count[2]  = { 1, 1 };
    hsize_t start[2]  = { samples_, 0 };
    hsize_t stride[2] = { 1, 1 };
    hsize_t block[2]  = { 1, 1 };
    ds_scalar_.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    hsize_t block2[2] = { 1, dimension };
    ds_vector_.selectHyperslab(H5S_SELECT_SET, count, start, stride, block2);

    hsize_t dim_vector[1] = { dimension };
    hsize_t dim_scalar[1] = { 1 };
    H5::DataSpace ds_vector(1, dim_vector);
    H5::DataSpace ds_scalar(1, dim_scalar);

    H5::DataType dt(H5::PredType::NATIVE_DOUBLE);

    // write thermodynamic equilibrium properties
    std::vector<H5::DataSet>::iterator dset = dset_.begin();
    (dset++)->write(&en_pot, dt, ds_scalar, ds_scalar_);
    (dset++)->write(&en_kin, dt, ds_scalar, ds_scalar_);
    (dset++)->write(&en_tot, dt, ds_scalar, ds_scalar_);
    (dset++)->write(&temp, dt, ds_scalar, ds_scalar_);
    (dset++)->write(&press, dt, ds_scalar, ds_scalar_);
    (dset++)->write(&v_cm, dt, ds_vector, ds_vector_);

    samples_++;
}

} // namespace mdsim

#endif /* ! MDSIM_ENERGY_HPP */
