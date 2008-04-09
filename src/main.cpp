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
#include "vector2d.hpp"
#include "vector3d.hpp"
#include "gsl_rng.hpp"
#include "ljfluid.hpp"
#include "mdsim.hpp"
using namespace std;


int main(int argc, char **argv)
{
    rng::gsl::gfsr4 rng;
#ifdef DIM_3D
    mdsim::ljfluid<vector3d<double> > fluid(100);
    mdsim::mdsim<vector3d<double> > sim;
#else
    mdsim::ljfluid<vector2d<double> > fluid(100);
    mdsim::mdsim<vector2d<double> > sim;
#endif

    rng.set(123);

    fluid.density(0.05);
    fluid.timestep(0.005);
    fluid.temperature(1., rng);

    for (size_t i = 1; i <= 10000; i++) {
	sim.step(fluid);

	if (i % 10) continue;

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
