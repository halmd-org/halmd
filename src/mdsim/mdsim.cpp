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

#include <boost/bind.hpp>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include "log.hpp"
#include "mdsim.hpp"
#include "signal.hpp"
#include "timer.hpp"

namespace mdsim
{

/**
 * initialize MD simulation program
 */
mdsim::mdsim(options const& opts) : opts(opts)
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
    // initialize random number generator with seed
    if (opts.rng_seed().empty()) {
	LOG("obtaining 32-bit integer seed from /dev/random");
	unsigned int seed;
	try {
	    std::ifstream rand;
	    rand.exceptions(std::ifstream::eofbit|std::ifstream::failbit|std::ifstream::badbit);
	    rand.open("/dev/random");
	    rand.read(reinterpret_cast<char*>(&seed), sizeof(seed));
	    rand.close();
	}
	catch (std::ifstream::failure const& e) {
	    throw std::logic_error(std::string("failed to read from /dev/random: ") + e.what());
	}
	fluid.rng(seed);
    }
    else {
	fluid.rng(opts.rng_seed().value());
    }

    if (!opts.trajectory_sample().empty()) {
	trajectory<false> traj;
	// open trajectory input file
	traj.open(opts.trajectory_input_file().value());
	// read trajectory sample and set system state
	fluid.restore(boost::bind(&trajectory<false>::read, boost::ref(traj), _1, _2, opts.trajectory_sample().value()));
	// close trajectory input file
	traj.close();
    }
    else {
	// arrange particles on a face-centered cubic (fcc) lattice
	fluid.lattice();
    }

    if (opts.trajectory_sample().empty() || !opts.temperature().defaulted()) {
	// set system temperature according to Maxwell-Boltzmann distribution
	fluid.temperature(opts.temperature().value());
    }
    // set simulation timestep
    fluid.timestep(opts.timestep().value());

    LOG("number of equilibration steps: " << opts.equilibration_steps().value());

#ifndef USE_BENCHMARK
    if (!opts.time().empty()) {
	// set total simulation time
	tcf.time(opts.time().value(), opts.timestep().value());
    }
    else {
	// set total number of simulation steps
	tcf.steps(opts.steps().value(), opts.timestep().value());
    }
    // set sample rate for lowest block level
    tcf.sample_rate(opts.sample_rate().value());
    // set block size
    tcf.block_size(opts.block_size().value());
    // set maximum number of samples per block
    tcf.max_samples(opts.max_samples().value());
    // set q-vectors for spatial Fourier transformation
    tcf.q_values(opts.q_values().value(), fluid.box());
#endif
}

/**
 * run MD simulation program
 */
void mdsim::operator()()
{
    // handle non-lethal POSIX signals to allow for a partial simulation run
    signal_handler signal;
    // measure elapsed realtime
    real_timer timer;

    // performance data
    prf.open(opts.output_file_prefix().value() + ".prf");
#ifndef USE_BENCHMARK
    prf.attrs() << fluid << tcf;
#else
    prf.attrs() << fluid;
#endif

    if (opts.equilibration_steps().value()) {
	LOG("starting equilibration");
	timer.start();
	for (iterator_timer<uint64_t> step = 0; step < opts.equilibration_steps().value(); ++step) {
	    // MD simulation step
	    fluid.mdstep();

	    // check whether a runtime estimate has finished
	    if (step.elapsed() > 0) {
		LOG("estimated remaining runtime: " << step);
		step.clear();
		// schedule next remaining runtime estimate
		step.set(TIME_ESTIMATE_INTERVAL);
	    }

	    // process signal event
	    if (*signal) {
		LOG_WARNING("trapped signal " << signal << " at simulation step " << *step);

		if (*signal == SIGUSR1) {
		    // schedule runtime estimate now
		    step.set(0);
		}
		else if (*signal == SIGINT || *signal == SIGTERM) {
		    LOG_WARNING("aborting equilibration");
		    signal.clear();
		    break;
		}
		signal.clear();
	    }
	}
	timer.stop();
	LOG("finished equilibration");
	LOG("total equilibration runtime: " << timer);
    }
    // sample performance counters
    prf.sample(fluid.times());
    // commit HDF5 performance datasets
    prf.commit();

#ifndef USE_BENCHMARK
    // time correlation functions
    tcf.open(opts.output_file_prefix().value() + ".tcf");
    tcf.attrs() << fluid << tcf;
    // trajectory file writer
    traj.open(opts.output_file_prefix().value() + ".trj", fluid.particles());
    traj.attrs() << fluid << tcf;
    // thermodynamic equilibrium properties
    tep.open(opts.output_file_prefix().value() + ".tep");
    tep.attrs() << fluid << tcf;

    // schedule first disk flush
    alarm(FLUSH_TO_DISK_INTERVAL);

    LOG("starting MD simulation");
    timer.start();

    for (iterator_timer<uint64_t> step = 0; step < tcf.steps(); ++step) {
	// check if sample is acquired for given simulation step
	if (tcf.sample(*step)) {
	    bool flush = false;
	    // simulation time
	    const double time = *step * fluid.timestep();
	    // sample time correlation functions
	    fluid.sample(boost::bind(&correlation::sample, boost::ref(tcf), _1, _2, *step, boost::ref(flush)));
	    // sample thermodynamic equilibrium properties
	    fluid.sample(boost::bind(&energy::sample, boost::ref(tep), _2, _3, _4, fluid.density(), time));
	    // sample trajectory
	    if (opts.dump_trajectories().value() || *step == 0) {
		fluid.sample(boost::bind(&trajectory<true>::sample, boost::ref(traj), _1, _2, time));
		if (*step == 0) {
		    traj.flush();
		}
	    }
	    // acquired maximum number of samples for a block level
	    if (flush) {
		// sample performance counters
		prf.sample(fluid.times());
		// write partial results to HDF5 files and flush to disk
		tcf.flush();
		if (opts.dump_trajectories().value())
		    traj.flush();
		tep.flush();
		prf.flush();
		LOG("flushed HDF5 buffers to disk");
		// schedule remaining runtime estimate
		step.clear();
		step.set(TIME_ESTIMATE_WAIT_AFTER_BLOCK);
		// schedule next disk flush
		alarm(FLUSH_TO_DISK_INTERVAL);
	    }
	}
	// MD simulation step
	fluid.mdstep();

	// check whether a runtime estimate has finished
	if (step.elapsed() > 0) {
	    LOG("estimated remaining runtime: " << step);
	    step.clear();
	    // schedule next remaining runtime estimate
	    step.set(TIME_ESTIMATE_INTERVAL);
	}

	// process signal event
	if (*signal) {
	    if (*signal != SIGALRM) {
		LOG_WARNING("trapped signal " << signal << " at simulation step " << *step);
	    }
	    if (*signal == SIGUSR1) {
		// schedule runtime estimate now
		step.set(0);
		signal.clear();
	    }
	    else if (*signal == SIGHUP || *signal == SIGALRM) {
		// sample performance counters
		prf.sample(fluid.times());
		// write partial results to HDF5 files and flush to disk
		tcf.flush();
		if (opts.dump_trajectories().value())
		    traj.flush();
		tep.flush();
		prf.flush();
		LOG("flushed HDF5 buffers to disk");
		// schedule next disk flush
		alarm(FLUSH_TO_DISK_INTERVAL);
	    }
	    else if (*signal == SIGINT || *signal == SIGTERM) {
		LOG_WARNING("aborting simulation");
		signal.clear();
		break;
	    }
	    signal.clear();
	}
    }

    // save last phase space sample
    fluid.sample(boost::bind(&trajectory<true>::sample, boost::ref(traj), _1, _2, tcf.steps() * fluid.timestep()));

    // sample performance counters
    prf.sample(fluid.times());
    // commit HDF5 performance datasets
    prf.commit();

    timer.stop();
    LOG("finished MD simulation");
    LOG("total MD simulation runtime: " << timer);

    // cancel previously scheduled disk flush
    alarm(0);
    // close HDF5 output files
    tcf.close();
    traj.close();
    tep.close();
#endif
    prf.close();
}

} // namespace mdsim
