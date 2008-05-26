/* Lennard-Jones fluid simulation with cell lists and Verlet neighbour lists
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

#include <algorithm>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include "H5param.hpp"
#include "exception.hpp"
#include "gsl_rng.hpp"
#include "log.hpp"
#include "options.hpp"


#define foreach BOOST_FOREACH

namespace mdsim
{

/**
 * MD simulation particle
 */
template <typename T>
struct particle
{
    /** particle position */
    T r;
    /** particle velocity */
    T v;
    /** particle number tag */
    unsigned int n;

    /** particle force */
    T f;
    /** particle neighbours list */
    std::vector<particle<T>* > neighbour;

    particle(T const& r, unsigned int n) : r(r), n(n) {}
    particle(T const& r, T const& v, unsigned int n) : r(r), v(v), n(n) {}
    particle() {}
};

/**
 * Lennard-Jones fluid simulation with cell lists and Verlet neighbour lists
 */
template <unsigned dimension, typename T>
class ljfluid
{
public:
    typedef typename std::list<particle<T> > cell_list;
    typedef typename cell_list::iterator cell_list_iterator;
    typedef typename cell_list::const_iterator cell_list_const_iterator;
    typedef typename boost::array<int, dimension> cell_index;

public:
    /** initialize fixed simulation parameters */
    ljfluid();
    /** set number of particles */
    void particles(unsigned int value);
    /** set system state from phase space sample */
    void particles(phase_space_point<std::vector<T> > const& state);
    /** set particle density */
    void density(double value);
    /** set periodic box length */
    void box(double value);
    /** initialize cell lists */
    void init_cell();
    /** set simulation timestep */
    void timestep(double value);

    /** initialize random number generator with seed */
    void rng(unsigned int seed);
    /** initialize random number generator from state */
    void rng(rng::gsl::gfsr4::state_type const& state);
    /** place particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(double value);

    /** returns number of particles */
    unsigned int const& particles() const;
    /** returns number of cells per dimension */
    unsigned int const& cells() const;
    /** returns particle density */
    double const& density() const;
    /** returns periodic box length */
    double const& box() const;
    /** returns cell length */
    double const& cell_length();
    /** returns simulation timestep */
    double const& timestep() const;

    /** copy ljfluid parameters to global simulation parameters */
    void copy_param(H5param& param) const;

    /** MD simulation step */
    void mdstep();
    /** sample trajectory */
    template <typename V> void sample(V& visitor) const;

private:
    /** update cell lists */
    void update_cells();
    /** returns cell list which a particle belongs to */
    cell_list& compute_cell(T const& r);
    /** update neighbour lists */
    void update_neighbours();
    /** update neighbour list of particle */
    void compute_cell_neighbours(particle<T>& p, cell_list& c);
    /** compute Lennard-Jones forces */
    void compute_forces();
    /** first leapfrog step of integration of equations of motion */
    void leapfrog_half();
    /** second leapfrog step of integration of equations of motion */
    void leapfrog_full();

private:
    /** number of particles */
    unsigned int npart;
    /** particle density */
    double density_;
    /** periodic box length */
    double box_;
    /** number of cells per dimension */
    unsigned int ncell;
    /** cell length */
    double cell_length_;
    /** simulation timestep */
    double timestep_;
    /** cutoff distance for shifted Lennard-Jones potential */
    double r_cut;
    /** neighbour list skin */
    double r_skin;
    /** cutoff distance with neighbour list skin */
    double r_cut_skin;

    /** cell lists */
    boost::multi_array<cell_list, dimension> cell;
    /** particles sorted by particle number */
    phase_space_point<std::vector<T> > part;
    /** potential energy per particle */
    double en_pot_;
    /** virial theorem force sum */
    double virial_;

    /** random number generator */
    rng::gsl::gfsr4 rng_;
    /** squared cutoff distance */
    double rr_cut;
    /** potential energy at cutoff distance */
    double en_cut;
    /** squared cutoff distance with neighbour list skin */
    double rr_cut_skin;
    /** sum over maximum velocity magnitudes since last neighbour lists update */
    double v_max_sum;
};

/**
 * initialize fixed simulation parameters
 */
template <unsigned dimension, typename T>
ljfluid<dimension, T>::ljfluid()
{
    // fixed cutoff distance for shifted Lennard-Jones potential
    r_cut = 2.5;
    LOG("potential cutoff distance: " << r_cut);

    // squared cutoff distance
    rr_cut = r_cut * r_cut;
    // potential energy at cutoff distance
    double rri_cut = 1. / rr_cut;
    double r6i_cut = rri_cut * rri_cut * rri_cut;
    en_cut = 4. * r6i_cut * (r6i_cut - 1.);

    // neighbour list skin
    r_skin = 0.3;
    LOG("neighbour list skin: " << r_skin);

    // cutoff distance with neighbour list skin
    r_cut_skin = r_skin + r_cut;
    // squared cutoff distance with neighbour list skin
    rr_cut_skin = r_cut_skin * r_cut_skin;
}

/**
 * set number of particles in system
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::particles(unsigned int value)
{
    if (value < 1) {
	throw exception("number of particles must be non-zero");
    }
    npart = value;
    LOG("number of particles: " << npart);

    try {
	part.r.resize(npart);
	part.R.resize(npart);
	part.v.resize(npart);
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate phase space state");
    }
}

/**
 * set system state from phase space sample
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::particles(phase_space_point<std::vector<T> > const& state)
{
    assert(state.r.size() == npart);
    assert(state.v.size() == npart);

    // copy particle positions to sorted particle list
    std::copy(state.r.begin(), state.r.end(), part.r.begin());
    std::copy(state.R.begin(), state.R.end(), part.R.begin());
    // copy particle velocities to sorted particle list
    std::copy(state.v.begin(), state.v.end(), part.v.begin());

    for (unsigned int i = 0; i < npart; ++i) {
	// add particle to appropriate cell list
	compute_cell(part.r[i]).push_back(particle<T>(part.r[i], part.v[i], i));
    }

    // update cell lists
    update_cells();
    // update Verlet neighbour lists
    update_neighbours();
    // reset sum over maximum velocity magnitudes to zero
    v_max_sum = 0.;
    // reconstruct forces for first leapfrog half step
    compute_forces();
}
     
/**
 * set particle density
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::density(double value)
{
    density_ = value;
    LOG("particle density: " << density_);

    // derive periodic box length
    box_ = std::pow(npart / density_, 1. / dimension);
    LOG("periodic box length: " << box_);
}

/**
 * set periodic box length
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::box(double value)
{
    box_ = value;
    LOG("periodic box length: " << box_);

    // derive particle density
    density_ = npart / std::pow(box_, 1. * dimension);
    LOG("particle density: " << density_);
}

/**
 * initialize cell lists
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::init_cell()
{
    // number of cells per dimension
    ncell = std::floor(box_ / r_cut_skin);
    LOG("number of cells per dimension: " << ncell);

    if (ncell < 3) {
	throw exception("requires at least 3 cells per dimension");
    }

    // derive cell length from integer number of cells per dimension
    cell_length_ = box_ / ncell;
    LOG("cell length: " << cell_length_);

    try {
#ifdef DIM_3D
	cell.resize(boost::extents[ncell][ncell][ncell]);
#else
	cell.resize(boost::extents[ncell][ncell]);
#endif
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate initial cell lists");
    }
}

/**
 * set simulation timestep
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::timestep(double value)
{
    timestep_ = value;
    LOG("simulation timestep: " << timestep_);
}

/**
 * initialize random number generator with seed
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::rng(unsigned int seed)
{
    rng_.set(seed);
    LOG("initializing random number generator with seed: " << seed);
}

/**
 * initialize random number generator from state
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::rng(rng::gsl::gfsr4::state_type const& state)
{
    rng_.restore(state);
    LOG("restoring random number generator from state");
}

/**
 * place particles on a face-centered cubic (fcc) lattice
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::lattice()
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
	// add particle to appropriate cell list
	compute_cell(r).push_back(particle<T>(r, i));
	// copy position to sorted particle list
	part.r[i] = r;
	part.R[i] = r;
    }
}

/**
 * set system temperature according to Maxwell-Boltzmann distribution
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::temperature(double value)
{
    // center of mass velocity
    T v_cm = 0.;
    // maximum squared velocity
    double vv_max = 0.;

    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle<T>& p, *it) {
	    // generate random Maxwell-Boltzmann distributed velocity
	    rng_.gaussian(p.v.x, p.v.y, value);
#ifdef DIM_3D
	    // Box-Muller transformation strictly generates 2 variates at once
	    rng_.gaussian(p.v.y, p.v.z, value);
#endif
	    v_cm += p.v;

	    // initialize force to zero for first leapfrog half step
	    p.f = 0.;
	}
    }

    v_cm /= npart;

    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle<T>& p, *it) {
	    // set center of mass velocity to zero
	    p.v -= v_cm;
	    // copy velocity to sorted particle list
	    part.v[p.n] = p.v;

	    vv_max = std::max(vv_max, p.v * p.v);
	}
    }

    // initialize sum over maximum velocity magnitudes since last neighbour lists update
    v_max_sum = std::sqrt(vv_max);
}

/**
 * returns number of particles
 */
template <unsigned dimension, typename T>
unsigned int const& ljfluid<dimension, T>::particles() const
{
    return npart;
}

/**
 * returns number of cells per dimension
 */
template <unsigned dimension, typename T>
unsigned int const& ljfluid<dimension, T>::cells() const
{
    return ncell;
}

/**
 * returns particle density
 */
template <unsigned dimension, typename T>
double const& ljfluid<dimension, T>::density() const
{
    return density_;
}

/**
 * returns periodic box length
 */
template <unsigned dimension, typename T>
double const& ljfluid<dimension, T>::box() const
{
    return box_;
}

/**
 * returns cell length
 */
template <unsigned dimension, typename T>
double const& ljfluid<dimension, T>::cell_length()
{
    return cell_length_;
}

/**
 * returns simulation timestep
 */
template <unsigned dimension, typename T>
double const& ljfluid<dimension, T>::timestep() const
{
    return timestep_;
}

/**
 * copy ljfluid parameters to global simulation parameters
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::copy_param(H5param& param) const
{
    // positional coordinate dimension
    param.dimension(dimension);
    // number of particles
    param.particles(npart);
    // number of cells per dimension
    param.cells(ncell);
    // particle density
    param.density(density_);
    // periodic box length
    param.box_length(box_);
    // cell length
    param.cell_length(cell_length_);
    // simulation timestep
    param.timestep(timestep_);
    // cutoff distance
    param.cutoff_distance(r_cut);
}

/**
 * update cell lists
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::update_cells()
{
    // FIXME particle may be inspected twice if moved ahead in cell order
    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	for (cell_list_iterator it2 = it->begin(); it2 != it->end(); ) {
	    // store old cell list iterator
	    cell_list_iterator p(it2);
	    // advance cell list iterator to next particle
	    it2++;
	    // compute cell which particle belongs to
	    cell_list& c = compute_cell(p->r);
	    // check whether cells are not identical
	    if (&c != &(*it)) {
		// move particle from old to new cell
		c.splice(c.end(), *it, p);
	    }
	}
    }
}

/**
 * returns cell list which a particle belongs to
 */
template <unsigned dimension, typename T>
typename ljfluid<dimension, T>::cell_list& ljfluid<dimension, T>::compute_cell(T const& r)
{
    T idx = trunc((r - floor(r / box_) * box_) / cell_length_);
#ifdef DIM_3D
    return cell[(int)(idx.x)][(int)(idx.y)][(int)(idx.z)];
#else
    return cell[(int)(idx.x)][(int)(idx.y)];
#endif
}

/**
 * update neighbour lists
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::update_neighbours()
{
#ifdef DIM_3D
    for (unsigned int x = 0; x < ncell; ++x) {
	for (unsigned int y = 0; y < ncell; ++y) {
	    for (unsigned int z = 0; z < ncell; ++z) {
		foreach (particle<T>& p, cell[x][y][z]) {
		    // empty neighbour list of particle
		    p.neighbour.clear();

		    boost::array<cell_index, 14> neighbour = {{
			{{  0,  0,  0 }},
			// visit half of neighbour cells (13 of 26) due to pair potential
			{{  0, -1,  0 }},
			{{ -1, -1,  0 }},
			{{ -1,  0,  0 }},
			{{ -1, +1,  0 }},
			{{  0, -1, -1 }},
			{{ -1, -1, -1 }},
			{{ -1,  0, -1 }},
			{{ -1, +1, -1 }},
			{{  0, -1, +1 }},
			{{ -1, -1, +1 }},
			{{ -1,  0, +1 }},
			{{ -1, +1, +1 }},
			{{  0,  0, -1 }},
		    }};

		    // update neighbour list of particle
		    foreach (cell_index const& idx, neighbour) {
			compute_cell_neighbours(p, cell[(x + ncell + idx[0]) % ncell][(y + ncell + idx[1]) % ncell][(z + ncell + idx[2]) % ncell]);
		    }
		}
	    }
	}
    }
#else
    for (unsigned int x = 0; x < ncell; ++x) {
	for (unsigned int y = 0; y < ncell; ++y) {
	    foreach (particle<T>& p, cell[x][y]) {
		// empty neighbour list of particle
		p.neighbour.clear();

		boost::array<cell_index, 5> neighbour = {{
		    {{  0,  0 }},
		    // visit half of neighbour cells (4 of 8) due to pair potential
		    {{  0, -1 }},
		    {{ -1, -1 }},
		    {{ -1,  0 }},
		    {{ -1, +1 }},
		}};

		// update neighbour list of particle
		foreach (cell_index const& idx, neighbour) {
		    compute_cell_neighbours(p, cell[(x + ncell + idx[0]) % ncell][(y + ncell + idx[1]) % ncell]);
		}
	    }
	}
    }
#endif
}

/**
 * update neighbour list of particle
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::compute_cell_neighbours(particle<T>& p, cell_list& c)
{
    foreach (particle<T>& pp, c) {
	// skip identical particles
	if (&p == &pp)
	    continue;

	// particle distance vector
	T r = p.r - pp.r;
	// enforce periodic boundary conditions
	r -= round(r / box_) * box_;
	// squared particle distance
	double rr = r * r;

	// enforce cutoff distance with neighbour list skin
	if (rr >= rr_cut_skin)
	    continue;

	// add particle to neighbour list
	p.neighbour.push_back(&pp);
    }
}

/**
 * compute Lennard-Jones forces
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::compute_forces()
{
    // initialize particle forces to zero
    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle<T>& p, *it) {
	    p.f = 0.;
	}
    }

    // potential energy
    en_pot_ = 0.;
    // virial equation sum
    virial_ = 0.;

    // iterate over all particles
    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle<T>& p, *it) {
	    // calculate pairwise Lennard-Jones force with neighbour particles
	    foreach (particle<T>* pp, p.neighbour) {
		// particle distance vector
		T r = p.r - pp->r;
		// enforce periodic boundary conditions
		r -= round(r / box_) * box_;
		// squared particle distance
		double rr = r * r;

		// enforce cutoff distance
		if (rr >= rr_cut)
		    continue;

		// compute Lennard-Jones force in reduced units
		double rri = 1. / rr;
		double r6i = rri * rri * rri;
		double fval = 48. * rri * r6i * (r6i - 0.5);

		// add force contribution to both particles
		p.f += r * fval;
		pp->f -= r * fval;

		// add contribution to potential energy
		en_pot_ += 4. * r6i * (r6i - 1.) - en_cut;
		// add contribution to virial equation sum
		virial_ += rr * fval;
	    }
	}
    }

    en_pot_ /= npart;
    virial_ /= npart;
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::leapfrog_half()
{
    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle<T>& p, *it) {
	    // half step velocity
	    p.v += p.f * (timestep_ / 2.);
	    // full step position
	    p.r += p.v * timestep_;
	    // copy position to sorted particle list
	    part.r[p.n] = p.r - floor(p.r / box_) * box_;
	    part.R[p.n] = p.r;
	}
    }
}

/**
 * second leapfrog step of integration of equations of motion
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::leapfrog_full()
{
    // maximum squared velocity
    double vv_max = 0.;

    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle<T>& p, *it) {
	    // full step velocity
	    p.v = p.v + p.f * (timestep_ / 2.);
	    // copy velocity to sorted particle list
	    part.v[p.n] = p.v;

	    vv_max = std::max(vv_max, p.v * p.v);
	}
    }

    // add to sum over maximum velocity magnitudes since last neighbour lists update
    v_max_sum += std::sqrt(vv_max);
}

/**
 * MD simulation step
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::mdstep()
{
    // calculate particle positions
    leapfrog_half();
    // update cell lists
    update_cells();

    if (v_max_sum * timestep_ > r_skin / 2.) {
	// update Verlet neighbour lists
	update_neighbours();
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0.;
    }

    // calculate Lennard-Jones forces
    compute_forces();
    // calculate velocities
    leapfrog_full();
}

/**
 * sample trajectory
 */
template <unsigned dimension, typename T>
template <typename V>
void ljfluid<dimension, T>::sample(V& visitor) const
{
    visitor.sample(part, en_pot_, virial_);
}

} // namespace mdsim

#undef foreach

#endif /* ! MDSIM_LJFLUID_HPP */
