/* Molecular Dynamics simulation of a Lennard-Jones fluid
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

#include <cuda_wrapper.hpp>
#include <stdint.h>
#include "autocorrelation.hpp"
#include "energy.hpp"
#include "ljfluid.hpp"
#include "options.hpp"
#include "rand48.hpp"
#include "trajectory.hpp"


namespace mdsim
{

/**
 * Molecular Dynamics simulation of a Lennard-Jones fluid
 */
template <unsigned dimension, typename T>
class mdsim
{
public:
    mdsim(options const& opts) : opts_(opts) {};
    void operator()();

private:
    /** program options */
    options const& opts_;
};


/**
 * run MD simulation
 */
template <unsigned dimension, typename T>
void mdsim<dimension, T>::operator()()
{
    //
    // initialize Lennard Jones fluid simulation
    //

    ljfluid<dimension, T> fluid;

    // set number of particles in system
    fluid.particles(opts_.particles());
    // set number of CUDA execution threads
    fluid.threads(opts_.threads());
    // initialize random number generator with seed
    fluid.rng(opts_.rngseed());

    // set particle density
    fluid.density(opts_.density());
    // arrange particles on a face-centered cubic (fcc) lattice
    fluid.lattice();
    // set system temperature according to Maxwell-Boltzmann distribution
    fluid.temperature(opts_.temperature());
    // set simulation timestep
    fluid.timestep(opts_.timestep());

    //
    // initialize trajectory sample visitors
    //

    // autocorrelation functions
    autocorrelation<dimension, T> tcf(opts_);
    // thermodynamic equilibrium properties
    energy<dimension, cuda::host::vector<T> > tep(opts_);
    // trajectory writer
    trajectory<dimension, cuda::host::vector<T> > traj(opts_);

    //
    // run MD simulation
    //

    // stream first MD simulation step on GPU
    fluid.mdstep();

    for (uint64_t step = 0; step < opts_.steps(); ++step) {
	// copy MD simulation step results from GPU to host
	fluid.sample();
	// stream next MD simulation step on GPU
	fluid.mdstep();

	// sample autocorrelation functions
	fluid.sample(tcf);
	// sample thermodynamic equilibrium properties
	fluid.sample(tep);
	// sample trajectory
	fluid.sample(traj);
    }

    //
    // write trajectory sample visitor results
    //

    // write autocorrelation function results to HDF5 file
    tcf.write(opts_.correlations_output_file(), fluid.timestep());
    // write thermodynamic equilibrium properties to HDF5 file
    tep.write(opts_.energy_output_file());
}

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
