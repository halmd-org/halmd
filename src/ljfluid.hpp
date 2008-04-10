/* Simulate a Lennard-Jones fluid with naive N-squared algorithm
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

#ifndef MDSIM_LJFLUID_HPP
#define MDSIM_LJFLUID_HPP

#include <vector>
#include <iostream>
#include <math.h>


namespace mdsim
{

/**
 * MD simulation particle
 */
template <typename T>
struct particle
{
    /** n-dimensional particle coordinates */
    T pos;
    /** n-dimensional particle velocity */
    T vel;
    /** n-dimensional force acting upon particle */
    T force;
};


/**
 * Simulate a Lennard-Jones fluid with naive N-squared algorithm
 */
template <typename T>
class ljfluid_base
{
public:
    typedef typename std::vector<particle<T> >::iterator iterator;
    typedef typename std::vector<particle<T> >::const_iterator const_iterator;

protected:
    /** number of particles in periodic box */
    size_t npart;
    /** particles */
    std::vector<particle<T> > part;

    /** particles per n-dimensional volume */
    double density_;
    /** MD simulation timestep */
    double timestep_;
    /** cutoff distance for shifted Lennard-Jones potential */
    double r_cut;

    /** periodic box length */
    double box;
    /** squared cutoff distance */
    double rr_cut;
    /** potential energy at cutoff distance */
    double en_cut;

public:
    /**
     * initialize Lennard-Jones fluid with given particle number
     */
    ljfluid_base(size_t npart) : npart(npart), part(npart)
    {
	// fixed cutoff distance for shifted Lennard-Jones potential
	// Frenkel
	r_cut = 2.5;
	// Rapaport
	//r_cut = pow(2., 1. / 6.);

	// squared cutoff distance
	rr_cut = r_cut * r_cut;

	// potential energy at cutoff distance
	double rri_cut = 1. / rr_cut;
	double r6i_cut = rri_cut * rri_cut * rri_cut;
	en_cut = 4. * r6i_cut * (r6i_cut - 1.);
    }

    /**
     * get number of particles in periodic box
     */
    size_t particles() const
    {
	return npart;
    }

    /**
     * get simulation timestep
     */
    double timestep()
    {
	return timestep_;
    }

    /**
     * set simulation timestep
     */
    void timestep(double val)
    {
	timestep_ = val;
    }

    /**
     * set temperature
     */
    template <typename rng_type>
    void temperature(double temp, rng_type& rng)
    {
	init_velocities(temp, rng);
	init_forces();
    }

    /**
     * MD simulation step
     */
    void step(double& en_pot, double& virial, T& vel_cm, double& vel2_sum)
    {
	leapfrog_half();
	compute_forces(en_pot, virial);
	leapfrog_full(vel_cm, vel2_sum);
    }

    /**
     * write particle coordinates and velocities to output stream
     */
    void trajectories(std::ostream& os) const
    {
	const_iterator it;
	size_t i;

	for (it = part.begin(), i = 1; it != part.end(); ++it, ++i) {
	    os << i << "\t" << it->pos << "\t" << it->vel << std::endl;
	}
	os << std::endl << std::endl;
    }

private:
    /*
     * first leapfrog half-step in integration of equations of motion
     */
    void leapfrog_half()
    {
	for (iterator it = part.begin(); it != part.end(); ++it) {
	    // half step velocity
	    it->vel += it->force * (timestep_ / 2.);
	    // full step coordinates
	    it->pos += it->vel * timestep_;
	    // enforce periodic boundary conditions
	    it->pos -= floor(it->pos / box) * box;
	}
    }

    /*
     * second leapfrog step in integration of equations of motion
     */
    void leapfrog_full(T& vel_cm, double& vel2_sum)
    {
	// center of mass velocity
	vel_cm = 0.;
	// squared velocities sum
	vel2_sum = 0.;

	for (iterator it = part.begin(); it != part.end(); ++it) {
	    // full step velocity
	    it->vel = it->vel + (timestep_ / 2) * it->force;
	    // velocity center of mass
	    vel_cm += it->vel;
	    // total kinetic energy
	    vel2_sum += it->vel * it->vel;
	}

	vel_cm /= npart;
	vel2_sum /= npart;
    }

    /*
     * compute pairwise Lennard-Jones forces
     */
    void compute_forces(double& en_pot, double& virial)
    {
	// potential energy
	en_pot = 0.;
	// virial equation sum
	virial = 0.;

	for (iterator it = part.begin(); it != part.end(); ++it) {
	    it->force = 0.;
	}

	for (iterator it = part.begin(); it != part.end(); ++it) {
	    for (iterator it2 = it; ++it2 != part.end(); ) {
		// particle distance vector
		T r = it->pos - it2->pos;
		// enforce periodic boundary conditions
		r -= round(r / box) * box;
		// squared particle distance
		double rr = r * r;

		// enforce cutoff distance
		if (rr >= rr_cut) continue;

		// compute Lennard-Jones force in reduced units
		double rri = 1. / rr;
		double r6i = rri * rri * rri;
		double fval = 48. * rri * r6i * (r6i - 0.5);

		it->force += r * fval;
		it2->force -= r * fval;

		en_pot += 4. * r6i * (r6i - 1.) - en_cut;
		virial += rr * fval;
	    }
	}

	en_pot /= npart;
	virial /= npart;
    }

    /**
     * generate random n-dimensional velocity vectors with uniform magnitude
     */
    template <typename rng_type>
    void init_velocities(double temp, rng_type& rng)
    {
	// velocity magnitude
	double vel_mag = sqrt(T::dim() * temp);
	// center of mass velocity
	T vel_cm = 0.;

	for (iterator it = part.begin(); it != part.end(); ++it) {
	    rng.unit_vector(it->vel);
	    it->vel *= vel_mag;
	    vel_cm += it->vel;
	}

	// set center of mass velocity to zero
	vel_cm /= npart;
	for (iterator it = part.begin(); it != part.end(); ++it) {
	    it->vel -= vel_cm;
	}
    }

    /**
     * set n-dimensional force vectors to zero
     */
    void init_forces()
    {
	for (iterator it = part.begin(); it != part.end(); ++it) {
	    it->force = 0.;
	}
    }
};


template <typename T>
class ljfluid;


/**
 * 2-dimensional Lennard-Jones fluid
 */
template <typename T>
class ljfluid<vector2d<T> > : public ljfluid_base<vector2d<T> >
{
private:
    typedef ljfluid_base<vector2d<T> > base_;

public:
    ljfluid(size_t npart) : base_(npart)
    {
    }

    /**
     * get particle density
     */
    double density() const
    {
	return base_::density_;
    }

    /**
     * set particle density
     */
    void density(double density_)
    {
	// particle density
	this->density_ = density_;
	// periodic box length
	this->box = sqrt(base_::npart / density_);

	init_lattice();
    }

private:
    /**
     * place particles on a 2-dimensional square lattice
     */
    void init_lattice()
    {
	typename base_::iterator it;
	size_t i;

	// number of particles along one lattice dimension
	size_t n = (size_t) ceil(sqrt(base_::npart));
	// lattice distance
	double a = base_::box / n;

	for (it = base_::part.begin(), i = 0; it != base_::part.end(); ++it, ++i) {
	    it->pos = vector2d<T>(i % n + 0.5, i / n + 0.5) * a;
	}
    }
};


/**
 * 3-dimensional Lennard-Jones fluid
 */
template <typename T>
class ljfluid<vector3d<T> > : public ljfluid_base<vector3d<T> >
{
private:
    typedef ljfluid_base<vector3d<T> > base_;

public:
    ljfluid(size_t npart) : ljfluid_base<vector3d<T> >(npart)
    {
    }

    /**
     * get particle density
     */
    double density() const
    {
	return base_::density_;
    }

    /**
     * set particle density
     */
    void density(double density_)
    {
	// particle density
	base_::density_ = density_;
	// periodic box length
	base_::box = cbrt(base_::npart / density_);

	init_lattice();
    }

private:
    /**
     * place particles on a 3-dimensional simple cubic lattice
     */
    void init_lattice()
    {
	typename base_::iterator it;
	size_t i;

	// number of particles along one lattice dimension
	size_t n = (size_t) ceil(cbrt(base_::npart));
	// lattice distance
	double a = base_::box / n;

	for (it = base_::part.begin(), i = 0; it != base_::part.end(); ++it, ++i) {
	    it->pos = vector3d<T>(i % n + 0.5, i / n % n + 0.5, i / n / n + 0.5) * a;
	}
    }
};

} // namespace mdsim

#endif /* ! MDSIM_LJFLUID_HPP */
