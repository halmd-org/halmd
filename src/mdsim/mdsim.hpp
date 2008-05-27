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
#include <boost/bind.hpp>
#include <stdint.h>
#include <vector>
#include "H5param.hpp"
#include "autocorrelation.hpp"
#include "energy.hpp"
#include "ljfluid.hpp"
#include "log.hpp"
#include "options.hpp"
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
    /** initialize MD simulation program */
    mdsim(options const& opts);
    /** run MD simulation program */
    void operator()();

private:
    /** program options */
    options const& opts;
    /** global simulation parameters */
    H5param param;
    /** Lennard-Jones fluid simulation */
    ljfluid<dimension, T> fluid;
};

/**
 * initialize MD simulation program
 */
template <unsigned dimension, typename T>
mdsim<dimension, T>::mdsim(options const& opts) : opts(opts)
{
    // initialize Lennard Jones fluid simulation
    if (!opts.trajectory_input_file().empty()) {
	// open trajectory input file
	trajectory<dimension, T, false> traj;
	traj.open(opts.trajectory_input_file().value());
	// read global simulation parameters
	traj.read(param);

	// set number of particles in system
	fluid.particles(param.particles());

	if (!opts.box_length().empty()) {
	    // set simulation box length
	    fluid.box(opts.box_length().value());
	}
	else {
	    // set particle density
	    fluid.density(opts.density().defaulted() ? param.density() : opts.density().value());
	}

	// initialize cell lists
	fluid.init_cell();
	// set simulation timestep
	fluid.timestep(opts.timestep().defaulted() ? param.timestep() : opts.timestep().value());

	// read trajectory sample and set system state
	fluid.state(boost::bind(&trajectory<dimension, T, false>::read, boost::ref(traj), _1, _2, opts.trajectory_sample().value()));
	// close trajectory input file
	traj.close();

	// initialize random number generator with seed
	fluid.rng(opts.rng_seed().value());

	if (!opts.temperature().empty()) {
	    // set system temperature according to Maxwell-Boltzmann distribution
	    fluid.temperature(opts.temperature().value());
	}
    }
    else {
	// set number of particles
	fluid.particles(opts.particles().value());
	// initialize random number generator with seed
	fluid.rng(opts.rng_seed().value());

	if (!opts.box_length().empty()) {
	    // set simulation box length
	    fluid.box(opts.box_length().value());
	}
	else {
	    // set particle density
	    fluid.density(opts.density().value());
	}

	// initialize cell lists
	fluid.init_cell();
	// set simulation timestep
	fluid.timestep(opts.timestep().value());

	// arrange particles on a face-centered cubic (fcc) lattice
	fluid.lattice();
	// set system temperature according to Maxwell-Boltzmann distribution
	fluid.temperature(opts.temperature().value());
    }
}

/**
 * run MD simulation program
 */
template <unsigned dimension, typename T>
void mdsim<dimension, T>::operator()()
{
    // autocorrelation functions
    autocorrelation<dimension, T> tcf(opts);
    // thermodynamic equilibrium properties
    energy<dimension, T> tep(opts);
    // trajectory writer
    trajectory<dimension, T> traj;

    // open HDF5 output files
    traj.open(opts);

    // collect global simulation parameters
    fluid.copy_param(param);
    tcf.copy_param(param);

    // write global simulation parameters to HDF5 output files
    tcf.write_param(param);
    tep.write_param(param);
    traj.write(param);

    // sample trajectory
    fluid.sample(boost::bind(&trajectory<dimension, T>::sample, boost::ref(traj), _1, _2, 0));

    for (uint64_t step = 0; step < opts.steps().value(); ++step) {
	// MD simulation step
	fluid.mdstep();

	// sample autocorrelation functions
	fluid.sample(boost::bind(&autocorrelation<dimension, T>::sample, boost::ref(tcf), _1, _2));
	// sample thermodynamic equilibrium properties
	fluid.sample(boost::bind(&energy<dimension, T>::sample, boost::ref(tep), _2, _3, _4));
	// sample trajectory
	fluid.sample(boost::bind(&trajectory<dimension, T>::sample, boost::ref(traj), _1, _2, step + 1));
    }

    // write autocorrelation function results to HDF5 file
    tcf.write(fluid.timestep());
    // write thermodynamic equilibrium properties to HDF5 file
    tep.write();
    // close HDF5 output files
    traj.close();
}

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
