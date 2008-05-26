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

#include <H5Cpp.h>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <cuda_wrapper.hpp>
#include <stdint.h>
#include "H5param.hpp"
#include "autocorrelation.hpp"
#include "energy.hpp"
#include "ljfluid.hpp"
#include "log.hpp"
#include "options.hpp"
#include "rand48.hpp"
#include "trajectory.hpp"
#include "version.h"


#define foreach BOOST_FOREACH

namespace mdsim
{

/**
 * Molecular Dynamics simulation of a Lennard-Jones fluid
 */
template <unsigned dimension, typename T>
class mdsim
{
public:
    typedef mdstep_sample<cuda::host::vector<T>, cuda::host::vector<float> > mdstep_sample_type;

public:
    mdsim(options const& opts) : opts_(opts) {};

    /** run MD simulation */
    void operator()();

private:
    /** program options */
    options const& opts_;
    /** Lennard-Jones fluid simulation */
    ljfluid<dimension, T> fluid;
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

    if (!opts_.trajectory_input_file().empty()) {
	// read trajectory sample
	phase_space_point<std::vector<T> > sample;
	trajectory<dimension, mdstep_sample_type>::read(opts_.trajectory_input_file().value(), opts_.trajectory_sample().value(), sample);
	// set system state from phase space sample
	fluid.particles(sample);

	// set number of CUDA execution threads
	fluid.threads(opts_.threads().value());
	// initialize random number generator with seed
	fluid.rng(opts_.rng_seed().value());

	if (!opts_.box_length().empty()) {
	    // set simulation box length
	    fluid.box(opts_.box_length().value());
	}
	else {
	    // set particle density
	    fluid.density(opts_.density().value());
	}

	if (!opts_.temperature().empty()) {
	    // set system temperature according to Maxwell-Boltzmann distribution
	    fluid.temperature(opts_.temperature().value());
	}
	// set simulation timestep
	fluid.timestep(opts_.timestep().value());
    }
    else {
	// set number of particles in system
	fluid.particles(opts_.particles().value());
	// set number of CUDA execution threads
	fluid.threads(opts_.threads().value());
	// initialize random number generator with seed
	fluid.rng(opts_.rng_seed().value());

	if (!opts_.box_length().empty()) {
	    // set simulation box length
	    fluid.box(opts_.box_length().value());
	}
	else {
	    // set particle density
	    fluid.density(opts_.density().value());
	}

	// arrange particles on a face-centered cubic (fcc) lattice
	fluid.lattice();
	// set system temperature according to Maxwell-Boltzmann distribution
	fluid.temperature(opts_.temperature().value());
	// set simulation timestep
	fluid.timestep(opts_.timestep().value());
    }

    //
    // initialize trajectory sample visitors
    //

    // autocorrelation functions
    autocorrelation<dimension, mdstep_sample_type> tcf(opts_);
    // thermodynamic equilibrium properties
    energy<dimension, mdstep_sample_type> tep(opts_);
    // trajectory writer
    trajectory<dimension, mdstep_sample_type> traj(opts_);

    // global simulation parameters
    H5param param;

    // gather simulation parameters
    fluid.copy_param(param);
    tcf.copy_param(param);
    // write simulation parameters to HDF5 output files
    tcf.write_param(param);
    tep.write_param(param);
    traj.write_param(param);

    //
    // run MD simulation
    //

    // stream first MD simulation step on GPU
    fluid.mdstep();

    for (uint64_t step = 0; step < tcf.steps(); ++step) {
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
    tcf.write(fluid.timestep());
    // write thermodynamic equilibrium properties to HDF5 file
    tep.write();
}

#undef foreach

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
