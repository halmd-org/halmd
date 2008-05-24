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
 * Molecular Dynamics simulation of a Lennard-Jones fluid
 */
template <unsigned dimension, typename T>
class mdsim
{
public:
    mdsim(options const& opts) : opts(opts) {}
    void operator()();

private:
    /** program options */
    options const& opts;
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

    // set number of particles
    fluid.particles(opts.particles().value());
    // initialize random number generator with seed
    fluid.rng(opts.rng_seed().value());

    if (!opts.box().empty()) {
	// set simulation box length
	fluid.box(opts.box().value());
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

    //
    // initialize trajectory sample visitors
    //

    // autocorrelation functions
    autocorrelation<dimension, T> tcf(opts);
    // thermodynamic equilibrium properties
    energy<dimension, std::vector<T> > tep(opts);
    // trajectory writer
    trajectory<dimension, std::vector<T> > traj(opts);

    //
    // run MD simulation
    //

    // sample trajectory
    fluid.sample(traj);

    for (uint64_t step = 0; step < opts.steps().value(); ++step) {
	// MD simulation step
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
    tcf.write(opts.output_file_prefix().value() + ".tcf", fluid.timestep());
    // write thermodynamic equilibrium properties to HDF5 file
    tep.write(opts.output_file_prefix().value() + ".tep");
}

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
