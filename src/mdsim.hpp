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
#include <stdint.h>
#include "H5param.hpp"
#include "autocorrelation.hpp"
#include "energy.hpp"
#include "ljfluid.hpp"
#include "log.hpp"
#include "options.hpp"
#include "rand48.hpp"
#include "trajectory.hpp"


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
    /** initialize MD simulation program */
    mdsim(options const& opts);
    /** run MD simulation program */
    void operator()();

private:
    /** program options */
    options const& opts;
    /** Lennard-Jones fluid simulation */
    ljfluid<dimension, T> fluid;
    /** trajectory file writer */
    trajectory<dimension, T> traj;
};

/**
 * initialize MD simulation program
 */
template <unsigned dimension, typename T>
mdsim<dimension, T>::mdsim(options const& opts) : opts(opts)
{
    // boost::bind
    using namespace boost;

    // initialize Lennard Jones fluid simulation
    if (!opts.trajectory_input_file().empty()) {
	// open trajectory input file
	trajectory<dimension, T, false> traj;
	traj.open(opts.trajectory_input_file().value());
	// read global simulation parameters
	H5param param;
	traj.read(param);
	// read trajectory sample and set system state
	fluid.state(bind(&trajectory<dimension, T, false>::read, ref(traj), _1, _2, opts.trajectory_sample().value()));
	// close trajectory input file
	traj.close();

	// set number of CUDA execution threads
	fluid.threads(opts.threads().defaulted() ? param.threads() : opts.threads().value());
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

	if (!opts.temperature().defaulted()) {
	    // set system temperature according to Maxwell-Boltzmann distribution
	    fluid.temperature(opts.temperature().value());
	}
	// set simulation timestep
	fluid.timestep(opts.timestep().defaulted() ? param.timestep() : opts.timestep().value());
    }
    else {
	// set number of particles in system
	fluid.particles(opts.particles().value());
	// set number of CUDA execution threads
	fluid.threads(opts.threads().value());
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

	// arrange particles on a face-centered cubic (fcc) lattice
	fluid.lattice();
	// set system temperature according to Maxwell-Boltzmann distribution
	fluid.temperature(opts.temperature().value());
	// set simulation timestep
	fluid.timestep(opts.timestep().value());
    }
}

/**
 * run MD simulation program
 */
template <unsigned dimension, typename T>
void mdsim<dimension, T>::operator()()
{
    // boost::bind
    using namespace boost;

    // autocorrelation functions
    autocorrelation<dimension, T> tcf(opts);
    // thermodynamic equilibrium properties
    energy<dimension, T> tep(opts);
    // open HDF5 output files
    traj.open(opts);

    H5param param;
    // collect global simulation parameters
    fluid.copy_param(param);
    tcf.copy_param(param);
    // write global simulation parameters to HDF5 output files
    tcf.write_param(param);
    tep.write_param(param);
    traj.write(param);

    // stream first MD simulation program step on GPU
    fluid.mdstep();

    for (uint64_t step = 0; step < tcf.steps(); ++step) {
	// copy MD simulation program step results from GPU to host
	fluid.sample();
	// stream next MD simulation program step on GPU
	fluid.mdstep();

	// sample autocorrelation functions
	fluid.sample(bind(&autocorrelation<dimension, T>::sample, ref(tcf), _1));
	// sample thermodynamic equilibrium properties
	fluid.sample(bind(&energy<dimension, T>::sample, ref(tep), _1));
	// sample trajectory
	fluid.sample(bind(&trajectory<dimension, T>::sample, ref(traj), _1, step));
    }

    // write autocorrelation function results to HDF5 file
    tcf.write(fluid.timestep());
    // write thermodynamic equilibrium properties to HDF5 file
    tep.write();
    // close HDF5 output files
    traj.close();
}

#undef foreach

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
