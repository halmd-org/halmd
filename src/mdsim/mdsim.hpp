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

#include <boost/bind.hpp>
#include <stdint.h>
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
    // set positional coordinates dimension
    param.dimension(dimension);
    LOG("positional coordinates dimension: " << dimension);

    // initialize Lennard Jones fluid simulation
    if (!opts.trajectory_input_file().empty()) {
	trajectory<dimension, T, false> traj;
	// open trajectory input file
	traj.open(opts.trajectory_input_file().value());
	// read global simulation parameters
	traj.read(param);

	// set number of particles
	fluid.particles(param.particles());

	if (!opts.box_length().empty()) {
	    // set simulation box length
	    fluid.box(opts.box_length().value());
	}
	else {
	    // set particle density
	    if (!opts.density().defaulted()) {
		param.density(opts.density().value());
	    }
	    fluid.density(param.density());
	}

	// initialize cell lists
	fluid.init_cell();

	// read trajectory sample and set system state
	fluid.state(boost::bind(&trajectory<dimension, T, false>::read, boost::ref(traj), _1, _2, opts.trajectory_sample().value()));
	// close trajectory input file
	traj.close();

	// initialize random number generator with seed
	fluid.rng(opts.rng_seed().value());

	if (!opts.temperature().defaulted()) {
	    LOG_WARNING("discarding velocities from trajectory file");
	    // set system temperature according to Maxwell-Boltzmann distribution
	    fluid.temperature(opts.temperature().value());
	}

	// override parameters with non-default option values
	if (!opts.timestep().defaulted()) {
	    param.timestep(opts.timestep().value());
	}
	if (!opts.block_size().defaulted()) {
	    param.block_size(opts.block_size().value());
	}
	if (!opts.max_samples().defaulted()) {
	    param.max_samples(opts.max_samples().value());
	}
	if (!opts.steps().defaulted()) {
	    param.steps(opts.steps().value());
	}
    }
    else {
	// set number of particles in system
	fluid.particles(opts.particles().value());

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

	// initialize random number generator with seed
	fluid.rng(opts.rng_seed().value());
	// arrange particles on a face-centered cubic (fcc) lattice
	fluid.lattice();
	// set system temperature according to Maxwell-Boltzmann distribution
	fluid.temperature(opts.temperature().value());

	// gather parameters from option values
	param.timestep(opts.timestep().value());
	param.block_size(opts.block_size().value());
	param.max_samples(opts.max_samples().value());
	param.steps(opts.steps().value());
    }

    if (!opts.time().empty()) {
	// override number of simulation steps with total simulation time
	param.time(opts.time().value());
	param.steps(std::ceil(param.time() / param.timestep()));
    }
    else {
	// derive total simulation time from number of simulation steps
	param.time(param.steps() * param.timestep());
    }
    LOG("total simulation time: " << param.time());

    // gather number of particles
    param.particles(fluid.particles());
    // gather number of cells per dimension
    param.cells(fluid.cells());
    // gather particle density
    param.density(fluid.density());
    // gather simulation box length
    param.box_length(fluid.box());
    // gather cell length
    param.cell_length(fluid.cell_length());
    // gather potential cutoff distance
    param.cutoff_distance(fluid.cutoff_distance());

    // set simulation timestep
    fluid.timestep(param.timestep());

    // adjust maximum number of samples to simulation steps limit
    param.max_samples(std::min(param.max_samples(), param.steps()));
    LOG("maximum number of samples: " << param.max_samples());
}

/**
 * run MD simulation program
 */
template <unsigned dimension, typename T>
void mdsim<dimension, T>::operator()()
{
    // autocorrelation functions
    autocorrelation<dimension, T> tcf(param);
    // gather block shift
    param.block_shift(tcf.block_shift());
    // gather block count
    param.block_count(tcf.block_count());

    // trajectory file writer
    trajectory<dimension, T> traj(param);
    // thermodynamic equilibrium properties
    energy<dimension, T> tep(param);

    if (opts.dry_run().value()) {
	// abort now that all simulation parameters have been handled
	return;
    }

    // create trajectory output file
    traj.open(opts.output_file_prefix().value() + ".trj");
    // create correlations output file
    tcf.open(opts.output_file_prefix().value() + ".tcf");
    // create thermodynamic equilibrium properties output file
    tep.open(opts.output_file_prefix().value() + ".tep");

    // write global simulation parameters to HDF5 output files
    traj << param;
    tcf << param;
    tep << param;

    // sample initial trajectory
    fluid.sample(boost::bind(&trajectory<dimension, T>::sample, boost::ref(traj), _1, _2));

    LOG("starting MD simulation");

    for (uint64_t step = 0; step < param.steps(); ++step) {
	// MD simulation step
	fluid.mdstep();

	// sample autocorrelation functions
	fluid.sample(boost::bind(&autocorrelation<dimension, T>::sample, boost::ref(tcf), _1, _2));
	// sample thermodynamic equilibrium properties
	fluid.sample(boost::bind(&energy<dimension, T>::sample, boost::ref(tep), _2, _3, _4));
	// sample trajectory
	fluid.sample(boost::bind(&trajectory<dimension, T>::sample, boost::ref(traj), _1, _2));
    }

    LOG("finished MD simulation");

    // write autocorrelation function results to HDF5 file
    tcf.write();
    // write thermodynamic equilibrium properties to HDF5 file
    tep.write();

    // close HDF5 output files
    tcf.close();
    tep.close();
    traj.close();
}

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
