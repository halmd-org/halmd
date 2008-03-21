/* Integrate equations of motion using two-step leapfrog method
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

#ifndef _LEAPFROG_HPP
#define _LEAPFROG_HPP

#include <vector>
using namespace std;


/**
 * Integrate equations of motion using two-step leapfrog method
 *
 * The algorithm is described in
 *
 * D. C. Rapaport, The Art of Molecular Dynamics Simulation, 2004, p. 19
 */
template <typename V, typename T>
class leapfrog
{
private:
    /** n-dimensional particle coordinates */
    vector<V>& r;
    /** n-dimensional particle velocities */
    vector<V>& v;
    /** n-dimensional particles forces */
    const vector<V>& f;

    /** time step */
    T dt;
    /** periodic box length */
    T box;

public:
    leapfrog(vector<V>& r, vector<V>& v, const vector<V>& f): r(r), v(v), f(f)
    {
	// FIXME reasonable default values depend on other parameters
	set_box_length(100.);
	set_timestep(0.01);
    }

    /**
     * Set periodic box length
     */
    void set_box_length(T box)
    {
	this->box = box;
    }

    /**
     * Set timestep
     */
    void set_timestep(T dt)
    {
	this->dt = dt;
    }

    /**
     * First leapfrog step of integration of equations of motion
     */
    void first()
    {
	size_t i;

	for (i = 0; i < r.size(); i++) {
	    // half step velocity
	    v[i] += f[i] * (dt / 2);
	    // full step coordinates
	    r[i] += v[i] * dt;
	    // enforce periodic boundary conditions
	    r[i] -= floor(r[i] / box) * box;
	}
    }

    /**
     * Second leapfrog step of integration of equations of motion
     */
    void second(V& v_cm, T& en_kin)
    {
	size_t i;

	v_cm = 0.;
	en_kin = 0.;

	for (i = 0; i < r.size(); i++) {
	    // full step velocity
	    v[i] = v[i] + f[i] * (dt / 2);
	    // velocity center of mass
	    v_cm += v[i];
	    // total kinetic energy
	    en_kin += v[i] * v[i];
	}

	v_cm /= r.size();
	en_kin *= 0.5 / r.size();
    }
};

#endif /* ! _LEAPFROG_HPP */
