/* Calculate particle forces using Lennard-Jones potential
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

#ifndef _FORCE_HPP
#define _FORCE_HPP

#include <vector>
using namespace std;


/**
 * Calculate particle forces using Lennard-Jones potential
 *
 * The algorithm is described in
 *
 * D. Frenkel & B. Smit, Understanding Molecular Simulation, 2002, p. 68
 */
template <typename V, typename T>
class ljforce
{
private:
    /** n-dimensional particle coordinates */
    const vector<V>& r;
    /** n-dimensional particles forces */
    vector<V>& f;

    /** periodic box length */
    T box;
    /** cutoff distance */
    T r_cut;
    /** squared cutoff distance */
    T r_cut_sq;
    /** cutoff energy at cutoff distance */
    T en_cut;

public:
    ljforce(const vector<V>& r, vector<V>& f) : r(r), f(f)
    {
	// FIXME reasonable default values depend on other parameters
	set_box_length(100.);
	set_cutoff_distance(2.5); // from Frenkel & Smit, p.98
    }

    /**
     * Set periodic box length
     */
    void set_box_length(T box)
    {
	this->box = box;
    }

    /**
     * Set cutoff distance
     */
    void set_cutoff_distance(T r_cut)
    {
	this->r_cut = r_cut;

	// squared cutoff distance
	this->r_cut_sq = r_cut * r_cut;

	// compute cutoff energy at cutoff distance
	r_cut = 1. / r_cut_sq;
	r_cut = r_cut * r_cut * r_cut;
	this->en_cut = 4. * r_cut * (r_cut - 1.);
    }

    /**
     * Calculate Lennard-Jones forces
     */
    void operator()(T& en_pot)
    {
	size_t i, j;
	V d;
	T d_sq, d_inv_sq, d_inv_6;

	// zero particle forces
	for (i = 0; i < r.size(); i++) {
	    f[i] = 0.;
	}
	// zero potential energy
	en_pot = 0.;

	// compute Lennard-Jones force for all particle pairs
	for (i = 0; i < r.size() - 1; i++) {
	    for (j = i + 1; j < r.size(); j++) {
		// particle distance vector
		d = r[i] - r[j];
		// enforce periodic boundary conditions
		d -= round(d / box) * box;
		// squared particle distance
		d_sq = d * d;

		// enforce cutoff distance
		if (d_sq >= r_cut_sq) continue;

		// compute Lennard-Jones force in reduced units
		d_inv_sq = 1. / d_sq;
		d_inv_6 = d_inv_sq * d_inv_sq * d_inv_sq;
		d *= 48. * d_inv_sq * d_inv_6 * (d_inv_6 - 0.5);
		// force contributes to both particles
		f[i] += d;
		f[j] -= d;

		// potential energy
		en_pot += 4. * d_inv_6 * (d_inv_6 - 1.) - en_cut;
	    }
	}

	en_pot /= r.size();
    }
};

#endif /* ! _FORCE_HPP */
