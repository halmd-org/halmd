/* Molecular Dynamics Simulation of a Lennard-Jones fluid
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

#include <iostream>
#include <stdint.h>
#include <vector>
#include "autocorrelation.hpp"
#include "energy.hpp"
#include "exception.hpp"
#include "ljfluid.hpp"
#include "options.hpp"
#include "time.hpp"
#include "trajectory.hpp"
#include "vector2d.hpp"
#include "vector3d.hpp"
#include "version.h"
using namespace std;


int main(int argc, char **argv)
{
    mdsim::options opts;

    try {
	opts.parse(argc, argv);
    }
    catch (mdsim::options::exit_exception const& e) {
	return e.status();
    }

#ifdef DIM_3D
    mdsim::ljfluid<3, vector3d<double> > fluid(opts);
    // thermodynamic equilibrium properties
    mdsim::energy<3, std::vector<vector3d<double> > > tep(opts);
    mdsim::trajectory<3, std::vector<vector3d<double> > > traj(opts);
    mdsim::autocorrelation<3, vector3d<double> > tcf(opts);
#else
    mdsim::ljfluid<2, vector2d<double> > fluid(opts);
    // thermodynamic equilibrium properties
    mdsim::energy<2, std::vector<vector2d<double> > > tep(opts);
    mdsim::trajectory<2, std::vector<vector2d<double> > > traj(opts);
    mdsim::autocorrelation<2, vector2d<double> > tcf(opts);
#endif

    if (opts.steps() < tcf.min_samples()) {
	throw mdsim::exception("less simulation steps than minimum required number of samples");
    }

    try {
	fluid.density(opts.density());
	fluid.timestep(opts.timestep());
	fluid.temperature(opts.temp());
    }
    catch (string const& e) {
	cerr << PROGRAM_NAME ": " << e << endl;
	return EXIT_FAILURE;
    }

    mdsim::timer timer;
    timer.start();

    for (uint64_t i = 1; i <= opts.steps(); i++) {
	fluid.mdstep();

	fluid.sample(tcf);
	fluid.sample(tep);
	fluid.sample(traj);
    }

    tcf.write(opts.correlations_output_file(), opts.timestep());
    tep.write(opts.energy_output_file());

    timer.stop();
    cerr << "Elapsed time: " << (timer.elapsed() * 1.E3) << "ms" << endl;

    return EXIT_SUCCESS;
}
