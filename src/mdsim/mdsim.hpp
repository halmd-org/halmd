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
#include "autocorrelation.hpp"
#include "perf.hpp"
#include "block.hpp"
#include "energy.hpp"
#include "ljfluid.hpp"
#include "log.hpp"
#include "options.hpp"
#include "signal.hpp"
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
    /** Lennard-Jones fluid simulation */
    ljfluid<dimension, T> fluid;
    /** block algorithm parameters */
    block_param<dimension, T> block;
};

/**
 * initialize MD simulation program
 */
template <unsigned dimension, typename T>
mdsim<dimension, T>::mdsim(options const& opts) : opts(opts)
{
    LOG("positional coordinates dimension: " << dimension);

    // set number of particles in system
    fluid.particles(opts.particles().value());
    // set simulation box length or particle density
    if (opts.density().defaulted() && !opts.box_length().empty())
	fluid.box(opts.box_length().value());
    else
	fluid.density(opts.density().value());
    // initialize cell lists
    fluid.init_cell();

    if (!opts.trajectory_input_file().empty()) {
	trajectory<dimension, T, false> traj;
	// open trajectory input file
	traj.open(opts.trajectory_input_file().value());
	// read trajectory sample and set system state
	fluid.restore(boost::bind(&trajectory<dimension, T, false>::read, boost::ref(traj), _1, _2, opts.trajectory_sample().value()));
	// close trajectory input file
	traj.close();
    }
    else {
	// arrange particles on a face-centered cubic (fcc) lattice
	fluid.lattice();
    }

    if (opts.trajectory_input_file().empty() || !opts.temperature().defaulted()) {
	// initialize random number generator with seed
	fluid.rng(opts.rng_seed().value());
	// set system temperature according to Maxwell-Boltzmann distribution
	fluid.temperature(opts.temperature().value());
    }
    // set simulation timestep
    fluid.timestep(opts.timestep().value());

#ifndef USE_BENCHMARK
    if (!opts.time().empty()) {
	// set total simulation time
	block.time(opts.time().value(), opts.timestep().value());
    }
    else {
	// set total number of simulation steps
	block.steps(opts.steps().value(), opts.timestep().value());
    }
    // set block size
    block.block_size(opts.block_size().value());
    // set maximum number of samples per block
    block.max_samples(opts.max_samples().value());
#endif
}

/**
 * run MD simulation program
 */
template <unsigned dimension, typename T>
void mdsim<dimension, T>::operator()()
{
#ifndef USE_BENCHMARK
    // time correlation functions
    autocorrelation<dimension, T> tcf(block, fluid.box(), opts.q_values().value());
    tcf.open(opts.output_file_prefix().value() + ".tcf");
    tcf.attrs() << fluid << block << tcf;
    // trajectory file writer
    trajectory<dimension, T> traj(block);
    if (opts.dump_trajectories().value()) {
	traj.open(opts.output_file_prefix().value() + ".trj", fluid.particles());
	traj.attrs() << fluid << block << tcf;
    }
    // thermodynamic equilibrium properties
    energy<dimension, T> tep(block);
    tep.open(opts.output_file_prefix().value() + ".tep");
    tep.attrs() << fluid << block << tcf;
#endif
    // performance data
    perf<dimension, T> prf;
    prf.open(opts.output_file_prefix().value() + ".prf");
#ifndef USE_BENCHMARK
    prf.attrs() << fluid << block << tcf;
#else
    prf.attrs() << fluid;
#endif

    // handle signals
    signal_handler signal;

    LOG("starting MD simulation");

    for (uint64_t step = 0; step < block.steps(); ++step) {
	// abort simulation on signal
	if (signal.get()) {
	    LOG_WARNING("caught signal at simulation step " << step);
	    break;
	}

#ifndef USE_BENCHMARK
	// sample time correlation functions
	fluid.sample(boost::bind(&autocorrelation<dimension, T>::sample, boost::ref(tcf), _1, _2));
	// sample thermodynamic equilibrium properties
	fluid.sample(boost::bind(&energy<dimension, T>::sample, boost::ref(tep), _2, _3, _4, fluid.density(), fluid.timestep()));
	if (opts.dump_trajectories().value()) {
	    // sample trajectory
	    fluid.sample(boost::bind(&trajectory<dimension, T>::sample, boost::ref(traj), _1, _2, fluid.particles(), fluid.timestep()));
	}
#endif

	// MD simulation step
	fluid.mdstep();
    }

    LOG("finished MD simulation");

#ifndef USE_BENCHMARK
    // write time correlation function results to HDF5 file
    tcf.write();
    tcf.close();
    // write thermodynamic equilibrium properties to HDF5 file
    tep.write();
    tep.close();
    // close HDF5 trajectory output file
    if (opts.dump_trajectories().value()) {
	traj.close();
    }
#endif
    // write performance data to HDF5 file
    prf.write(fluid.times());
    prf.close();
}

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
