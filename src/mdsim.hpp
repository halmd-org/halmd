/* Molecular Dynamics simulation
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

#ifndef MDSIM_MDSIM_HPP
#define MDSIM_MDSIM_HPP

#include <stdint.h>
#include <vector>
#include "autocorrelation.hpp"
#include "energy.hpp"
#include "ljfluid.hpp"
#include "options.hpp"
#include "trajectory.hpp"


namespace mdsim
{

/**
 * Molecular Dynamics simulation
 */
template <unsigned dimension, typename T>
class mdsim
{
public:
    mdsim(options const& options);
    void run();

private:
    // Lennard Jones fluid simulation
    ljfluid<dimension, T> fluid_;
    // autocorrelation functions
    autocorrelation<dimension, T> tcf_;
    // thermodynamic equilibrium properties
    energy<dimension, std::vector<T> > tep_;
    // trajectory writer
    trajectory<dimension, std::vector<T> > traj_;

    // program options
    options const& opts_;
};


template <unsigned dimension, typename T>
mdsim<dimension, T>::mdsim(options const& opts) : fluid_(opts), tcf_(opts), tep_(opts), traj_(opts), opts_(opts)
{
    fluid_.density(opts_.density().value());
    fluid_.timestep(opts_.timestep().value());
    fluid_.temperature(opts_.temperature().value());
}


template <unsigned dimension, typename T>
void mdsim<dimension, T>::run()
{
    for (uint64_t step = 0; step < opts_.steps().value(); step++) {
	// MD simulation step
	fluid_.mdstep();

	// sample autocorrelation functions
	fluid_.sample(tcf_);
	// sample thermodynamic equilibrium properties
	fluid_.sample(tep_);
	// sample trajectory
	fluid_.sample(traj_);
    }

    // write autocorrelation function results to HDF5 file
    tcf_.write(opts_.output_file_prefix().value() + ".tcf", fluid_.timestep());
    // write thermodynamic equilibrium properties to HDF5 file
    tep_.write(opts_.output_file_prefix().value() + ".tep");
}

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
