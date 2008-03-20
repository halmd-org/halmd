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
#include <vector>
using namespace std;

#include "vector2d.hpp"
#include "gsl_rng.hpp"
#include "leapfrog.hpp"
#include "force.hpp"


int main(int argc, char **argv)
{
    unsigned int i, j;
    rng::gsl::gfsr4 rng;

    vector<vector2d<double> > r(100);
    vector<vector2d<double> > v(r.size());
    vector<vector2d<double> > f(r.size());

    vector2d<double> v_cm;
    double en_pot;
    double en_kin;

    ljforce<vector2d<double>, double> force(r, f);
    leapfrog<vector2d<double>, double> inteq(r, v, f);

    const double box = 100.;
    inteq.set_timestep(1E-2);
    //const double box = 50.;
    //inteq.set_timestep(1E-4);
    //const double box = 10.;
    //inteq.set_timestep(1E-9);

    force.set_box_length(box);
    inteq.set_box_length(box);

    // FIXME /dev/random
    rng.set(1);

    for (i = 0; i < r.size(); i++) {
	v[i].x = rng.get_uniform() - 0.5;
	v[i].y = rng.get_uniform() - 0.5;

	r[i].x = rng.get_uniform() * box;
	r[i].y = rng.get_uniform() * box;
    }

    // output gnuplot header
    cout << "# i\tr[i].x\tr[i].y\tf[i].x\tf[i].y" << endl;

    for (i = 0; i < 10000; i++) {
	// first leapfrog step of integration of equations of motion
	inteq.first();
	// force calculation
	force(en_pot);
	// second leapfrog step of integration of equations of motion
	inteq.second(v_cm, en_kin);

	if (i % 10) continue;
#if 1
	cout << "# en_pot(" << en_pot << ")" << endl;
	cout << "# en_kin(" << en_kin << ")" << endl;
	cout << "# en_tot(" << en_pot + en_kin << ")" << endl;
	cout << "# v_cm(" << v_cm.x << ", " << v_cm.y << ")" << endl;

	for (j = 0; j < f.size(); j++) {
	    cout << j << "\t" << r[j].x << "\t" << r[j].y << "\t";
	    cout << v[j].x << "\t" << v[j].y << "\t";
	    cout << f[j].x << "\t" << f[j].y << endl;
	}
	cout << endl << endl;
#endif
    }

    return EXIT_SUCCESS;
}
