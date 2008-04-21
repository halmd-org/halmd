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
#include "vector2d.hpp"
#include "vector3d.hpp"
#include "gsl_rng.hpp"
#include "ljfluid.hpp"
#include "mdsim.hpp"
#include "options.hpp"
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

    rng::gsl::gfsr4 rng;
#ifdef DIM_3D
    mdsim::ljfluid<vector3d<double> > fluid(opts.npart());
    mdsim::mdsim<vector3d<double> > sim;
#else
    mdsim::ljfluid<vector2d<double> > fluid(opts.npart());
    mdsim::mdsim<vector2d<double> > sim;
#endif

    rng.set(opts.rngseed());

    try {
	fluid.density(opts.density());
	fluid.timestep(opts.timestep());
	fluid.temperature(opts.temp(), rng);
    }
    catch (string const& e) {
	cerr << PROGRAM_NAME ": " << e << endl;
	return EXIT_FAILURE;
    }

    for (uint64_t i = 1; i <= opts.steps(); i++) {
	sim.step(fluid);

	if (i % opts.avgsteps())
	    continue;

	cout << "## steps(" << i << ")" << endl;

	cout << "# en_pot(" << sim.en_pot().mean() << ")" << endl;
	cout << "# sigma_en_pot(" << sim.en_pot().std() << ")" << endl;
	cout << "# en_kin(" << sim.en_kin().mean() << ")" << endl;
	cout << "# sigma_en_kin(" << sim.en_kin().std() << ")" << endl;
	cout << "# en_tot(" << sim.en_tot().mean() << ")" << endl;
	cout << "# sigma_en_tot(" << sim.en_tot().std() << ")" << endl;
	cout << "# temp(" << sim.temp().mean() << ")" << endl;
	cout << "# sigma_temp(" << sim.temp().std() << ")" << endl;
	cout << "# pressure(" << sim.pressure().mean() << ")" << endl;
	cout << "# sigma_pressure(" << sim.pressure().std() << ")" << endl;
	cout << "# vel_cm(" << sim.vel_cm().mean() << ")" << endl;
	cout << "# sigma_vel_cm(" << sim.vel_cm().std() << ")" << endl;

	fluid.trajectories(cout);

	sim.clear();
    }

    return EXIT_SUCCESS;
}
