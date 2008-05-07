/* Simulate a Lennard-Jones fluid with cell list algorithm
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

#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include "cell_array.hpp"
#include "trajectory.hpp"


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
    /** particle number */
    unsigned int tag;

    std::vector<particle<T>*> neighbour;

    particle() {}
    particle(T const& pos, unsigned int const& tag) : pos(pos), tag(tag) {}
};


/**
 * Simulate a Lennard-Jones fluid with naive N-squared algorithm
 */
template <unsigned int NDIM, typename T>
class ljfluid
{
protected:
    typedef typename std::list<particle<T> > list_type;
    typedef typename std::list<particle<T> >::iterator list_iterator;
    typedef typename std::list<particle<T> >::const_iterator list_const_iterator;

    typedef typename cell_array<list_type, T>::iterator cell_iterator;
    typedef typename cell_array<list_type, T>::const_iterator cell_const_iterator;

public:
    ljfluid(size_t npart);

    size_t particles() const;
    double timestep();
    void timestep(double val);
    double density() const;
    void density(double density_);
    double box() const;
    template <typename rng_type>
    void temperature(double temp, rng_type& rng);

    void step(double& en_pot, double& virial, T& vel_cm, double& vel2_sum);
    void trajectories(std::ostream& os) const;
    template <typename Y>
    void trajectories(trajectory<NDIM, Y>& traj) const;

private:
    void leapfrog_half();
    void leapfrog_full(T& vel_cm, double& vel2_sum);
    void update_cells();

    void compute_neighbours();
    void compute_cell_neighbours(particle<T>& p, list_type& cell);
    void compute_neighbour(particle<T>& p1, particle<T>& p2);

    void compute_forces(double& en_pot, double& virial);
    void compute_force(particle<T>& p1, particle<T>& p2, double& en_pot, double& virial);

    void init_cells();
    void lattice();
    template <typename rng_type>
    void init_velocities(double temp, rng_type& rng);
    void init_forces();

private:
    /** number of particles in periodic box */
    size_t npart;
    /** cell lists */
    cell_array<list_type, T> cells;
    /** number of cells along 1 dimension */
    size_t ncell;
    /** cell edge length */
    double cell_len;

    /** particles sorted by particle number */
    std::vector<T> part_;
    std::vector<T> vel_;
    std::vector<T> force_;

    /** particle density */
    double density_;
    /** periodic box length */
    double box_;
    /** MD simulation timestep */
    double timestep_;
    /** cutoff distance for shifted Lennard-Jones potential */
    double r_cut;

    /** squared cutoff distance */
    double rr_cut;
    /** potential energy at cutoff distance */
    double en_cut;
    
    /** neighbour list radius */
    double r_skin;
    double r_cut_skin;
    double rr_cut_skin;
    double vel_max_sum;
};


/**
 * initialize Lennard-Jones fluid with given particle number
 */
template <unsigned int NDIM, typename T>
ljfluid<NDIM, T>::ljfluid(size_t npart) : npart(npart), part_(npart), vel_(npart), force_(npart)
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

    r_skin = 0.3; // FIXME should depend on system size
    r_cut_skin = r_cut + r_skin;
    rr_cut_skin = r_cut_skin * r_cut_skin;
}

/**
 * get number of particles in periodic box
 */
template <unsigned int NDIM, typename T>
size_t ljfluid<NDIM, T>::particles() const
{
    return npart;
}

/**
 * get simulation timestep
 */
template <unsigned int NDIM, typename T>
double ljfluid<NDIM, T>::timestep()
{
    return timestep_;
}

/**
 * set simulation timestep
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::timestep(double val)
{
    timestep_ = val;
}

/**
 * get particle density
 */
template <unsigned int NDIM, typename T>
double ljfluid<NDIM, T>::density() const
{
    return density_;
}

/**
 * set particle density
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::density(double density_)
{
    // particle density
    this->density_ = density_;
    // periodic box length
    this->box_ = pow(npart / density_, 1.0 / T::dim());

    init_cells();
    lattice();
}

/**
 * get periodic box length
 */
template <unsigned int NDIM, typename T>
double ljfluid<NDIM, T>::box() const
{
    return box_;
}

/**
 * set temperature
 */
template <unsigned int NDIM, typename T>
template <typename rng_type>
void ljfluid<NDIM, T>::temperature(double temp, rng_type& rng)
{
    init_velocities(temp, rng);
    init_forces();
}

/**
 * MD simulation step
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::step(double& en_pot, double& virial, T& vel_cm, double& vel2_sum)
{
    // calculate coordinates
    leapfrog_half();
    // update cell lists
    update_cells();

    // update Verlet neighbour lists if necessary
    if ((2. * timestep_ * vel_max_sum) > r_skin) {
	compute_neighbours();
	vel_max_sum = 0.;
    }

    // calculate forces
    compute_forces(en_pot, virial);
    // calculate velocities
    leapfrog_full(vel_cm, vel2_sum);

    // sort particles by particle number for trajectory sampling
    for (cell_const_iterator cell = cells.begin(); cell != cells.end(); ++cell) {
	for (list_const_iterator it = cell->begin(); it != cell->end(); ++it) {
	    part_[it->tag] = it->pos;
	    vel_[it->tag] = it->vel;
	    force_[it->tag] = it->force;
	}
    }
}

/**
 * write particle coordinates and velocities to output stream
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::trajectories(std::ostream& os) const
{
    for (unsigned int i = 0; i < npart; ++i) {
	os << part_[i] << "\t" << vel_[i] << "\n";
    }
    os << std::endl << std::endl;
}

/**
 * write particle coordinates and velocities to binary HDF5 file
 */
template <unsigned int NDIM, typename T>
template <typename Y>
void ljfluid<NDIM, T>::trajectories(trajectory<NDIM, Y>& traj) const
{
    traj.write(part_, vel_, force_);
}

/**
 * first leapfrog half-step in integration of equations of motion
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::leapfrog_half()
{
    for (cell_iterator cell = cells.begin(); cell != cells.end(); ++cell) {
	for (list_iterator it = cell->begin(); it != cell->end(); ++it) {
	    // half step velocity
	    it->vel += it->force * (timestep_ / 2.);
	    // full step coordinates
	    it->pos += it->vel * timestep_;
	    // enforce periodic boundary conditions
	    it->pos -= floor(it->pos / box_) * box_;

	    // if particles are moving more than cutoff length in timestep,
	    // something is really fishy...
	    assert(((it->vel * it->vel) * timestep_ * timestep_) < rr_cut);
	}
    }
}

/**
 * second leapfrog step in integration of equations of motion
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::leapfrog_full(T& vel_cm, double& vel2_sum)
{
    // center of mass velocity
    vel_cm = 0.;
    // squared velocities sum
    vel2_sum = 0.;
    // maximum squared velocity
    double vel2_max = 0.;

    for (cell_iterator cell = cells.begin(); cell != cells.end(); ++cell) {
	for (list_iterator it = cell->begin(); it != cell->end(); ++it) {
	    // full step velocity
	    it->vel = it->vel + (timestep_ / 2) * it->force;
	    // velocity center of mass
	    vel_cm += it->vel;
	    // total kinetic energy
	    vel2_sum += it->vel * it->vel;
	    // maximum squared velocity
	    vel2_max = std::max(it->vel * it->vel, vel2_max);
	}
    }

    vel_cm /= npart;
    vel2_sum /= npart;
    vel_max_sum += sqrt(vel2_max);
}

/**
 * update cell lists
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::update_cells()
{
    for (cell_iterator cell = cells.begin(); cell != cells.end(); ++cell) {
	// FIXME particles may be visited twice if moved ahead in sequence
	for (list_iterator it = cell->begin(); it != cell->end(); ) {
	    list_iterator it_old = it;
	    ++it;

	    // update cell lists
	    list_type& cell_ = cells(it_old->pos / cell_len);
	    if (&cell_ != &(*cell)) {
		cell_.splice(cell_.end(), *cell, it_old);
	    }
	}
    }
}

/**
 * compute pairwise Lennard-Jones forces
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::compute_forces(double& en_pot, double& virial)
{
    // potential energy
    en_pot = 0.;
    // virial equation sum
    virial = 0.;

    for (cell_iterator cell = cells.begin(); cell != cells.end(); ++cell) {
	for (list_iterator it = cell->begin(); it != cell->end(); ++it) {
	    it->force = 0.;
	}
    }

    for (cell_iterator cell = cells.begin(); cell != cells.end(); ++cell) {
	for (list_iterator it = cell->begin(); it != cell->end(); ++it) {
	    for (typename std::vector<particle<T>*>::iterator it2 = it->neighbour.begin(); it2 != it->neighbour.end(); it2++) {
		compute_force(*it, **it2, en_pot, virial);
	    }
	}
    }

    en_pot /= npart;
    virial /= npart;
}

/**
 * FIXME
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::compute_neighbours()
{
    for (size_t i = 0; i < cells.size(); ++i) {
	for (size_t j = 0; j < cells[i].size(); ++j) {
#ifdef DIM_3D
	    for (size_t k = 0; k < cells[i][j].size(); ++k) {
		// FIXME This is royal prefix/postfix fun...
		for (list_iterator it = cells[i][j][k].begin(); it != cells[i][j][k].end(); ++it) {
		    // empty neighbour list
		    it->neighbour.clear();

		    // FIXME This is royal prefix/postfix fun...
		    for (list_iterator it2 = it; ++it2 != cells[i][j][k].end(); )
			compute_neighbour(*it, *it2);

		    // only half of neighbour cells need to be considered due to pair potential
		    compute_cell_neighbours(*it, cells[i][(j + ncell - 1) % ncell][k]);
		    compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][(j + ncell - 1) % ncell][k]);
		    compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][j][k]);
		    compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][(j + 1) % ncell][k]);

		    compute_cell_neighbours(*it, cells[i][(j + ncell - 1) % ncell][(k + ncell - 1) % ncell]);
		    compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][(j + ncell - 1) % ncell][(k + ncell - 1) % ncell]);
		    compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][j][(k + ncell - 1) % ncell]);
		    compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][(j + 1) % ncell][(k + ncell - 1) % ncell]);

		    compute_cell_neighbours(*it, cells[i][(j + ncell - 1) % ncell][(k + 1) % ncell]);
		    compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][(j + ncell - 1) % ncell][(k + 1) % ncell]);
		    compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][j][(k + 1) % ncell]);
		    compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][(j + 1) % ncell][(k + 1) % ncell]);

		    compute_cell_neighbours(*it, cells[i][j][(k + ncell - 1) % ncell]);
		}
	    }
#else
	    // FIXME This is royal prefix/postfix fun...
	    for (list_iterator it = cells[i][j].begin(); it != cells[i][j].end(); ++it) {
		// empty neighbour list
		it->neighbour.clear();

		// FIXME This is royal prefix/postfix fun...
		for (list_iterator it2 = it; ++it2 != cells[i][j].end(); )
		    compute_neighbour(*it, *it2);

		// only half of neighbour cells need to be considered due to pair potential
		compute_cell_neighbours(*it, cells[i][(j + ncell - 1) % ncell]);
		compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][(j + ncell - 1) % ncell]);
		compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][j]);
		compute_cell_neighbours(*it, cells[(i + ncell - 1) % ncell][(j + 1) % ncell]);
	    }
#endif
	}
    }
}

/**
 * FIXME
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::compute_cell_neighbours(particle<T>& p, list_type& cell)
{
    for (list_iterator it = cell.begin(); it != cell.end(); ++it) {
	compute_neighbour(p, *it);
    }
}

/**
 * FIXME
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::compute_neighbour(particle<T>& p1, particle<T>& p2)
{
    T r = p1.pos - p2.pos;
    // enforce periodic boundary conditions
    r -= round(r / box_) * box_;
    // squared particle distance
    double rr = r * r;

    // enforce cutoff distance
    if (rr >= rr_cut_skin) return;

    p1.neighbour.push_back(&p2);
}

/**
 * compute pairwise Lennard-Jones force
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::compute_force(particle<T>& p1, particle<T>& p2, double& en_pot, double& virial)
{
    T r = p1.pos - p2.pos;
    // enforce periodic boundary conditions
    r -= round(r / box_) * box_;
    // squared particle distance
    double rr = r * r;

    // enforce cutoff distance
    if (rr >= rr_cut) return;

    // compute Lennard-Jones force in reduced units
    double rri = 1. / rr;
    double r6i = rri * rri * rri;
    double fval = 48. * rri * r6i * (r6i - 0.5);

    p1.force += r * fval;
    p2.force -= r * fval;

    en_pot += 4. * r6i * (r6i - 1.) - en_cut;
    virial += rr * fval;
}

/**
 * initialize cell lists
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::init_cells()
{
    // number of cells per dimension
    ncell = size_t(floor(box_ / r_cut_skin));
    if (ncell < 3) {
	throw std::string("requires at least 3 cells per dimension");
    }
    // cell edge length (must be greater or equal to cutoff length)
    cell_len = box_ / ncell;

    cells.resize(ncell);
}

/**
 * place particles on a face centered cubic (FCC) lattice
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::lattice()
{
#ifdef DIM_3D
    // number of particles along 1 lattice dimension
    const unsigned int n = ceil(cbrt(npart / 4.));
#else
    // number of particles along 1 lattice dimension
    const unsigned int n = ceil(sqrt(npart / 2.));
#endif
    // lattice distance
    double a = box_ / n;

    for (unsigned int i = 0; i < npart; ++i) {
	T r(a);
#ifdef DIM_3D
	// compose primitive vectors from 1-dimensional index
	r.x *= ((i >> 2) % n) + ((i ^ (i >> 1)) & 1) / 2.;
	r.y *= ((i >> 2) / n % n) + (i & 1) / 2.;
	r.z *= ((i >> 2) / n / n) + (i & 2) / 4.;
#else
	// compose primitive vectors from 1-dimensional index
	r.x *= ((i >> 1) % n) + (i & 1) / 2.;
	r.y *= ((i >> 1) / n) + (i & 1) / 2.;
#endif
	cells(r / cell_len).push_back(particle<T>(r, i));
    }
}

/**
 * generate random n-dimensional Maxwell-Boltzmann distributed velocities
 */
template <unsigned int NDIM, typename T>
template <typename rng_type>
void ljfluid<NDIM, T>::init_velocities(double temp, rng_type& rng)
{
    // center of mass velocity
    T vel_cm = 0.;
    // maximum squared velocity
    double vel2_max = 0.;

    for (cell_iterator cell = cells.begin(); cell != cells.end(); ++cell) {
	for (list_iterator it = cell->begin(); it != cell->end(); ++it) {
	    rng.gaussian(it->vel.x, it->vel.y, temp);
#ifdef DIM_3D
	    // Box-Muller transformation strictly generates 2 variates at once
	    rng.gaussian(it->vel.y, it->vel.z, temp);
#endif
	    vel_cm += it->vel;
	}
    }

    // set center of mass velocity to zero
    vel_cm /= npart;
    for (cell_iterator cell = cells.begin(); cell != cells.end(); ++cell) {
	for (list_iterator it = cell->begin(); it != cell->end(); ++it) {
	    it->vel -= vel_cm;
	    vel2_max = std::max(it->vel * it->vel, vel2_max);
	}
    }

    vel_max_sum = sqrt(vel2_max);
}

/**
 * set n-dimensional force vectors to zero
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::init_forces()
{
    for (cell_iterator cell = cells.begin(); cell != cells.end(); ++cell) {
	for (list_iterator it = cell->begin(); it != cell->end(); ++it) {
	    it->force = 0.;
	}
    }
}

} // namespace mdsim

#endif /* ! MDSIM_LJFLUID_HPP */
