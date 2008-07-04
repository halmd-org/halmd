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

#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <cuda_wrapper.hpp>
#include <iomanip>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "correlation.hpp"
#include "energy.hpp"
#include "ljfluid.hpp"
#include "log.hpp"
#include "options.hpp"
#include "perf.hpp"
#include "rand48.hpp"
#include "signal.hpp"
#include "trajectory.hpp"


#define foreach BOOST_FOREACH

namespace mdsim
{

/**
 * Molecular Dynamics simulation of a Lennard-Jones fluid
 */
template <unsigned dimension, typename T, typename U>
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
    ljfluid<dimension, T, U> fluid;
#ifndef USE_BENCHMARK
    /** block correlations */
    correlation<dimension, T, U> tcf;
    /**  trajectory file writer */
    trajectory<dimension, T, U> traj;
    /** thermodynamic equilibrium properties */
    energy<dimension, T, U> tep;
#endif
    /** performance data */
    perf<dimension, T, U> prf;
};

/**
 * initialize MD simulation program
 */
template <unsigned dimension, typename T, typename U>
mdsim<dimension, T, U>::mdsim(options const& opts) : opts(opts)
{
    LOG("positional coordinates dimension: " << dimension);

    // set number of particles in system
    fluid.particles(opts.particles().value());
    // set simulation box length or particle density
    if (opts.density().defaulted() && !opts.box_length().empty())
	fluid.box(opts.box_length().value());
    else
	fluid.density(opts.density().value());
#ifdef USE_CELL
    // compute cell parameters
    fluid.cell_occupancy(opts.cell_occupancy().value());
#endif
    // set number of CUDA execution threads
    fluid.threads(opts.threads().value());

    if (!opts.trajectory_sample().empty()) {
	trajectory<dimension, T, U, false> traj;
	// open trajectory input file
	traj.open(opts.trajectory_input_file().value());
	// read trajectory sample and restore system state
	fluid.restore(boost::bind(&trajectory<dimension, T, U, false>::read, boost::ref(traj), _1, _2, opts.trajectory_sample().value()));
	// close trajectory input file
	traj.close();
    }
    else {
	// arrange particles on a face-centered cubic (fcc) lattice
	fluid.lattice();
    }

    if (opts.trajectory_sample().empty() || !opts.temperature().defaulted()) {
	// initialize random number generator with seed
	fluid.rng(opts.rng_seed().value());
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
template <unsigned dimension, typename T, typename U>
void mdsim<dimension, T, U>::operator()()
{
#ifndef USE_BENCHMARK
    // time correlation functions
    tcf.open(opts.output_file_prefix().value() + ".tcf");
    tcf.attrs() << fluid << tcf;
    // trajectory file writer
    if (opts.dump_trajectories().value()) {
	traj.open(opts.output_file_prefix().value() + ".trj", fluid.particles());
	traj.attrs() << fluid << tcf;
    }
    // thermodynamic equilibrium properties
    tep.open(opts.output_file_prefix().value() + ".tep");
    tep.attrs() << fluid << tcf;
#endif
    // performance data
    prf.open(opts.output_file_prefix().value() + ".prf");
#ifndef USE_BENCHMARK
    prf.attrs() << fluid << tcf;
#else
    prf.attrs() << fluid;
#endif

    if (opts.equilibration_steps().value()) {
	LOG("starting equilibration");
	for (uint64_t step = 0; step <= opts.equilibration_steps().value(); ++step) {
	    // copy previous MD simulation state from GPU to host
	    fluid.sample();
	    // stream next MD simulation program step on GPU
	    fluid.mdstep();
	    // synchronize MD simulation program step on GPU
	    fluid.synchronize();
	}
	LOG("finished equilibration");
    }

#ifndef USE_BENCHMARK
    // handle non-lethal POSIX signals to allow for a partial simulation run
    signal_handler signal;

    // elapsed real time measurement
    std::pair<timeval, uint64_t> tv0, tv1, tv2;
    gettimeofday(&tv0.first, NULL);
    tv0.second = 0;
    tv1 = tv0;
    // schedule initial estimate of remaining MD simulation runtime
    alarm(300);

    LOG("starting MD simulation");
    for (uint64_t step = 0; step <= tcf.steps(); ++step) {
	// copy previous MD simulation state from GPU to host
	fluid.sample();
	// stream next MD simulation program step on GPU
	fluid.mdstep();
	// check if sample is acquired for given simulation step
	if (tcf.sample(step)) {
	    // simulation time
	    const float time = step * fluid.timestep();
	    // sample time correlation functions
	    fluid.sample(boost::bind(&correlation<dimension, T, U>::sample, boost::ref(tcf), _2, _3, step));
	    // sample thermodynamic equilibrium properties
	    fluid.sample(boost::bind(&energy<dimension, T, U>::sample, boost::ref(tep), _3, _4, _5, fluid.density(), time));
	    // sample trajectory
	    if (opts.dump_trajectories().value()) {
		fluid.sample(boost::bind(&trajectory<dimension, T, U>::sample, boost::ref(traj), _1, _2, _3, time));
	    }
	}
	// synchronize MD simulation program step on GPU
	fluid.synchronize();

	if (signal.get() == SIGALRM || signal.get() == SIGUSR1) {
	    // estimate remaining MD simulation runtime
	    gettimeofday(&tv2.first, NULL);
	    timersub(&tv2.first, &tv1.first, &tv1.first);
	    tv2.second = step + 1;
	    tv1.second = tv2.second - tv1.second;
	    const float trem = (tv1.first.tv_sec + tv1.first.tv_usec * 1.E-6) * (tcf.steps() + 1 - step) / tv1.second;
	    if (trem < 60)
		LOG("estimated remaining runtime: " << std::fixed << std::setprecision(1) << trem << " s");
	    else if (trem < 3600)
		LOG("estimated remaining runtime: " << std::fixed << std::setprecision(1) << (trem / 60) << " min");
	    else
		LOG("estimated remaining runtime: " << std::fixed << std::setprecision(1) << (trem / 3600) << " h");
	    tv1 = tv2;
	    signal.clear();
	    // schedule next estimate
	    alarm(1800);
	}
	else if (signal.get()) {
	    // abort simulation on signal
	    LOG_WARNING("caught signal at simulation step " << step);
	    signal.clear();
	    break;
	}
    }
    LOG("finished MD simulation");

    // output total MD simulation runtime
    gettimeofday(&tv1.first, NULL);
    timersub(&tv1.first, &tv0.first, &tv1.first);
    const float ttot = tv1.first.tv_sec + tv1.first.tv_usec * 1.E-6;
    if (ttot < 60)
	LOG("total MD simulation runtime: " << std::fixed << std::setprecision(1) << ttot << " s");
    else if (ttot < 3600)
	LOG("total MD simulation runtime: " << std::fixed << std::setprecision(1) << (ttot / 60) << " min");
    else
	LOG("total MD simulation runtime: " << std::fixed << std::setprecision(1) << (ttot / 3600) << " h");

    // write time correlation function results to HDF5 file
    tcf.write();
    tcf.close();
    // close HDF5 trajectory output file
    if (opts.dump_trajectories().value()) {
	traj.close();
    }
    // write thermodynamic equilibrium properties to HDF5 file
    tep.write();
    tep.close();
#endif
    // write performance data to HDF5 file (includes equilibration phase)
    prf.write(fluid.times());
    prf.close();
}

#undef foreach

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
