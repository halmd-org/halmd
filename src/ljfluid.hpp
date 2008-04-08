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
#include <math.h>
#include "gsl_rng.hpp"


namespace mdsim
{

template <typename V>
class measure;


/**
 * Simulate a Lennard-Jones fluid with naive N-squared algorithm
 */
template <typename V>
class _ljfluid_base
{
    friend class measure<V>;

protected:
    /** number of particles in periodic box */
    size_t N;
    /** n-dimensional particle coordinates */
    std::vector<V> pos;
    /** n-dimensional particle velocities */
    std::vector<V> vel;
    /** n-dimensional particles forces */
    std::vector<V> force;

    /** random number generator */
    rng::gsl::rand48 rng;

    /** potential energy accumulator */
    std::vector<double> en_pot;
    /** kinetic energy accumulator */
    std::vector<double> en_kin;
    /** virial equation sum accumulator */
    std::vector<double> virial;
    /** center of mass velocity accumulator */
    std::vector<V> vel_cm;

    /** particles per n-dimensional volume */
    double rho;
    /** MD simulation timestep */
    double h;
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
    _ljfluid_base(size_t N) : N(N), pos(N), vel(N), force(N)
    {
	// fixed cutoff distance for shifted Lennard-Jones potential
	this->r_cut = 2.5;

	// squared cutoff distance
	this->rr_cut = r_cut * r_cut;

	// potential energy at cutoff distance
	double rri_cut = 1. / rr_cut;
	double r6i_cut = rri_cut * rri_cut * rri_cut;
	this->en_cut = 4. * r6i_cut * (r6i_cut - 1.);

	// FIXME
	rng.set(123);
    }

    /**
     * get number of particles in periodic box
     */
    size_t particles() const
    {
	return N;
    }

    /**
     * get particle density
     */
    double density() const
    {
	return rho;
    }

    /**
     * get simulation timestep
     */
    double timestep()
    {
	return h;
    }

    /**
     * set simulation timestep
     */
    void timestep(double h)
    {
	this->h = h;
    }

    /**
     * set temperature
     */
    void temperature(double temp)
    {
	init_velocities(temp);
	init_forces();
    }

    /**
     * MD simulation step
     */
    void step()
    {
	leapfrog_half();
	compute_forces();
	leapfrog_full();
    }

    /**
     * write particle coordinates and velocities to output stream
     */
    void trajectories(std::ostream& os)
    {
	for (size_t i = 0; i < N; i++) {
	    os << i + 1 << "\t" << pos[i] << "\t" << vel[i] << endl;
	}
	os << endl << endl;
    }

private:
    /*
     * first leapfrog half-step in integration of equations of motion
     */
    void leapfrog_half()
    {
	for (size_t i = 0; i < N; i++) {
	    // half step velocity
	    vel[i] += force[i] * (h / 2.);
	    // full step coordinates
	    pos[i] += vel[i] * h;
	    // enforce periodic boundary conditions
	    pos[i] -= floor(pos[i] / box) * box;
	}
    }

    /*
     * second leapfrog step in integration of equations of motion
     */
    void leapfrog_full()
    {
	// center of mass velocity
	V vel_cm = 0.;
	// kinetic energy
	double en_kin = 0.;

	for (size_t i = 0; i < N; i++) {
	    // full step velocity
	    vel[i] = vel[i] + (h / 2) * force[i];
	    // velocity center of mass
	    vel_cm += vel[i];
	    // total kinetic energy
	    en_kin += vel[i] * vel[i];
	}

	// accumulate center of mass velocity
	this->vel_cm.push_back(vel_cm / N);
	// accumulate kinetic energy
	this->en_kin.push_back(0.5 * en_kin / N);
    }

    /*
     * compute pairwise Lennard-Jones forces
     */
    void compute_forces()
    {
	// potential energy
	double en_pot = 0.;
	// virial equation sum
	double virial = 0.;

	for (size_t i = 0; i < N; i++) {
	    force[i] = 0.;
	}

	for (size_t i = 0; i < N - 1; i++) {
	    for (size_t j = i + 1; j < N; j++) {
		// particle distance vector
		V r = pos[i] - pos[j];
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

		force[i] += r * fval;
		force[j] -= r * fval;

		en_pot += 4. * r6i * (r6i - 1.) - en_cut;
		virial += rr * fval;
	    }
	}

	// accumulate potential energy
	this->en_pot.push_back(en_pot / N);
	// accumulate virial equation sum
	this->virial.push_back(virial / N);
    }

    /**
     * generate random n-dimensional velocity vectors with uniform magnitude
     */
    void init_velocities(double temp)
    {
	// velocity magnitude
	double vel_mag = sqrt(V::dim() * temp);
	// center of mass velocity
	V vel_cm = 0.;

	for (size_t i = 0; i < N; i++) {
	    rng.unit_vector(vel[i]);
	    vel[i] *= vel_mag;
	    vel_cm += vel[i];
	}

	// set center of mass velocity to zero
	vel_cm /= N;
	for (size_t i = 0; i < N; i++) {
	    vel[i] -= vel_cm;
	}
    }

    /**
     * set n-dimensional force vectors to zero
     */
    void init_forces()
    {
	for (size_t i = 0; i < N; i++) {
	    force[i] = 0.;
	}
    }
};


template <typename V>
class ljfluid;


/**
 * 2-dimensional Lennard-Jones fluid
 */
template <>
class ljfluid<vector2d<double> > : public _ljfluid_base<vector2d<double> >
{
private:
    typedef vector2d<double> V;

public:
    ljfluid(size_t N) : _ljfluid_base<V>(N)
    {
    }

    /**
     * set particle density
     */
    void density(double rho)
    {
	// particle density
	this->rho = rho;
	// periodic box length
	this->box = sqrt(this->N / rho);

	init_lattice();
    }

private:
    /**
     * place particles on a 2-dimensional square lattice
     */
    void init_lattice()
    {
	// number of particles along one lattice dimension
	size_t n = (size_t) ceil(sqrt(this->N));
	// lattice distance
	double a = this->box / n;

	for (size_t i = 0; i < this->N; i++) {
	    this->pos[i] = V(i % n + 0.5, i / n + 0.5) * a;
	}
    }
};


/**
 * 3-dimensional Lennard-Jones fluid
 */
template <>
class ljfluid<vector3d<double> > : public _ljfluid_base<vector3d<double> >
{
private:
    typedef vector3d<double> V;

public:
    ljfluid(size_t N) : _ljfluid_base<V>(N)
    {
    }

    /**
     * set particle density
     */
    void density(double rho)
    {
	// particle density
	this->rho = rho;
	// periodic box length
	this->box = cbrt(this->N / rho);

	init_lattice();
    }

private:
    /**
     * place particles on a 3-dimensional simple cubic lattice
     */
    void init_lattice()
    {
	// number of particles along one lattice dimension
	size_t n = (size_t) ceil(cbrt(this->N));
	// lattice distance
	double a = this->box / n;

	for (size_t i = 0; i < this->N; i++) {
	    this->pos[i] = V(i % n + 0.5, i / n % n + 0.5, i / n / n + 0.5) * a;
	}
    }
};

} // namespace mdsim

#endif /* ! MDSIM_LJFLUID_HPP */
