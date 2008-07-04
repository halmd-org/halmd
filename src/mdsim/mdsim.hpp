/* Event-based Molecular Dynamics simulation of hard spheres
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
#include <iomanip>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "correlation.hpp"
#include "energy.hpp"
#include "hardspheres.hpp"
#include "log.hpp"
#include "options.hpp"
#include "perf.hpp"
#include "signal.hpp"
#include "trajectory.hpp"


namespace mdsim
{

/**
 * Event-based Molecular Dynamics simulation of hard spheres
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
    /** hard spheres simulation */
    hardspheres<dimension, T> fluid;
#ifndef USE_BENCHMARK
    /** block correlations */
    correlation<dimension, T> tcf;
    /**  trajectory file writer */
    trajectory<dimension, T> traj;
    /** thermodynamic equilibrium properties */
    energy<dimension, T> tep;
#endif
    /** performance data */
    perf<dimension, T> prf;
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
    // set pair separation at which particle collision occurs
    fluid.pair_separation(opts.pair_separation().value());
    // set simulation box length or particle density
    if (opts.density().defaulted() && !opts.box_length().empty())
	fluid.box(opts.box_length().value());
    else
	fluid.density(opts.density().value());
    // initialize cells
    fluid.init_cell();

    if (!opts.trajectory_sample().empty()) {
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
    if (opts.trajectory_sample().empty() || !opts.temperature().defaulted()) {
	// initialize random number generator with seed
	fluid.rng(opts.rng_seed().value());
	// set system temperature according to Maxwell-Boltzmann distribution
	fluid.temperature(opts.temperature().value());
    }
    // initialize event list
    fluid.init_event_list();

    LOG("number of equilibration steps: " << opts.equilibration_steps().value());

#ifndef USE_BENCHMARK
    if (!opts.time().empty()) {
	// set total sample time
	tcf.time(opts.time().value(), opts.timestep().value());
    }
    else {
	// set total number of sample steps
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
template <unsigned dimension, typename T>
void mdsim<dimension, T>::operator()()
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
#ifndef USE_BENCHMARK
	// advance phase space state to given sample time
	fluid.mdstep(opts.equilibration_steps().value() * tcf.timestep());
#else
	for (uint64_t step = 0; step < opts.equilibration_steps().value(); ++step) {
	    // advance phase space state to given sample time
	    fluid.mdstep(step * opts.timestep().value());
	}
#endif
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
    alarm(60);

    LOG("starting MD simulation");
    for (uint64_t step = 0; step <= tcf.steps(); ++step) {
	// check if sample is acquired for given simulation step
	if (tcf.sample(step)) {
	    // sample time
	    const double time = step * tcf.timestep();
	    // sample time correlation functions
	    fluid.sample(boost::bind(&correlation<dimension, T>::sample, boost::ref(tcf), _2, _3, step));
	    // sample thermodynamic equilibrium properties
	    fluid.sample(boost::bind(&energy<dimension, T>::sample, boost::ref(tep), _3, _4, fluid.density(), tcf.timestep(), time));
	    // sample trajectory
	    if (opts.dump_trajectories().value()) {
		fluid.sample(boost::bind(&trajectory<dimension, T>::sample, boost::ref(traj), _1, _3, time));
	    }
	    // advance phase space state to given sample time
	    fluid.mdstep(step * tcf.timestep());
	}

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

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
