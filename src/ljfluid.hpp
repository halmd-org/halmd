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
 * Simulate a Lennard-Jones fluid with naive N-squared algorithm
 */
template <typename V>
class ljfluid_base
{
protected:
    /** number of particles in periodic box */
    size_t N;
    /** n-dimensional particle coordinates */
    std::vector<V> pos;
    /** n-dimensional particle velocities */
    std::vector<V> vel;
    /** n-dimensional particles forces */
    std::vector<V> force;

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
    ljfluid_base(size_t N) : N(N), pos(N), vel(N), force(N)
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
	return N;
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
    template <typename T>
    void temperature(double temp, T& rng)
    {
	init_velocities(temp, rng);
	init_forces();
    }

    /**
     * MD simulation step
     */
    void step(double& en_pot, double& virial, V& vel_cm, double& vel2_sum)
    {
	leapfrog_half();
	compute_forces(en_pot, virial);
	leapfrog_full(vel_cm, vel2_sum);
    }

    /**
     * write particle coordinates and velocities to output stream
     */
    void trajectories(std::ostream& os)
    {
	for (size_t i = 0; i < N; i++) {
	    os << i + 1 << "\t" << pos[i] << "\t" << vel[i] << std::endl;
	}
	os << std::endl << std::endl;
    }

private:
    /*
     * first leapfrog half-step in integration of equations of motion
     */
    void leapfrog_half()
    {
	for (size_t i = 0; i < N; i++) {
	    // half step velocity
	    vel[i] += force[i] * (timestep_ / 2.);
	    // full step coordinates
	    pos[i] += vel[i] * timestep_;
	    // enforce periodic boundary conditions
	    pos[i] -= floor(pos[i] / box) * box;
	}
    }

    /*
     * second leapfrog step in integration of equations of motion
     */
    void leapfrog_full(V& vel_cm, double& vel2_sum)
    {
	// center of mass velocity
	vel_cm = 0.;
	// squared velocities sum
	vel2_sum = 0.;

	for (size_t i = 0; i < N; i++) {
	    // full step velocity
	    vel[i] = vel[i] + (timestep_ / 2) * force[i];
	    // velocity center of mass
	    vel_cm += vel[i];
	    // total kinetic energy
	    vel2_sum += vel[i] * vel[i];
	}

	vel_cm /= N;
	vel2_sum /= N;
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

	en_pot /= N;
	virial /= N;
    }

    /**
     * generate random n-dimensional velocity vectors with uniform magnitude
     */
    template <typename T>
    void init_velocities(double temp, T& rng)
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
class ljfluid<vector2d<double> > : public ljfluid_base<vector2d<double> >
{
private:
    typedef vector2d<double> V;

public:
    ljfluid(size_t N) : ljfluid_base<V>(N)
    {
    }

    /**
     * get particle density
     */
    double density() const
    {
	return density_;
    }

    /**
     * set particle density
     */
    void density(double density_)
    {
	// particle density
	this->density_ = density_;
	// periodic box length
	this->box = sqrt(this->N / density_);

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
class ljfluid<vector3d<double> > : public ljfluid_base<vector3d<double> >
{
private:
    typedef vector3d<double> V;

public:
    ljfluid(size_t N) : ljfluid_base<V>(N)
    {
    }

    /**
     * get particle density
     */
    double density() const
    {
	return density_;
    }

    /**
     * set particle density
     */
    void density(double density_)
    {
	// particle density
	this->density_ = density_;
	// periodic box length
	this->box = cbrt(this->N / density_);

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
