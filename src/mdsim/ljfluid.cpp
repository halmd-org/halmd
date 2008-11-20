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

#include <algorithm>
#include <boost/foreach.hpp>
#include <cmath>
#include <iostream>
#include <sys/times.h>
#include "H5xx.hpp"
#include "config.hpp"
#include "exception.hpp"
#include "gsl_rng.hpp"
#include "ljfluid.hpp"
#include "log.hpp"
#include "perf.hpp"
#include "timer.hpp"

#define foreach BOOST_FOREACH

namespace mdsim
{

/**
 * set potential cutoff radius
 */
void ljfluid::cutoff_radius(double value)
{
    r_cut = value;
    LOG("potential cutoff radius: " << r_cut);

    // squared cutoff radius
    rr_cut = r_cut * r_cut;
    // potential energy at cutoff radius
    double rri_cut = 1. / rr_cut;
    double r6i_cut = rri_cut * rri_cut * rri_cut;
    en_cut = 4. * r6i_cut * (r6i_cut - 1.);

    // neighbour list skin
    r_skin = 0.5;
    LOG("neighbour list skin: " << r_skin);

    // cutoff radius with neighbour list skin
    r_cut_skin = r_skin + r_cut;
    // squared cutoff radius with neighbour list skin
    rr_cut_skin = r_cut_skin * r_cut_skin;
}

/**
 * set number of particles in system
 */
void ljfluid::particles(unsigned int value)
{
    if (value < 1) {
	throw exception("number of particles must be non-zero");
    }
    npart = value;
    LOG("number of particles: " << npart);

    try {
	part.r.resize(npart);
	part.v.resize(npart);
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate phase space state");
    }
}

/**
 * set system state from phase space sample
 */
void ljfluid::restore(trajectory_sample_visitor visitor)
{
    // set system state from phase space sample
    visitor(part.r, part.v);

    for (unsigned int i = 0; i < npart; ++i) {
	// add particle to appropriate cell list
	compute_cell(part.r[i]).push_back(particle(part.r[i], part.v[i], i));
    }

    // update Verlet neighbour lists
    update_neighbours();
    // reset sum over maximum velocity magnitudes to zero
    v_max_sum = 0.;
    // calculate forces, potential energy and virial equation sum
    compute_forces();
}

/**
 * set particle density
 */
void ljfluid::density(double value)
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
void ljfluid::box(double value)
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
void ljfluid::init_cell()
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
void ljfluid::timestep(double value)
{
    timestep_ = value;
    LOG("simulation timestep: " << timestep_);
}

/**
 * initialize random number generator with seed
 */
void ljfluid::rng(unsigned int seed)
{
    rng_.set(seed);
    LOG("initializing random number generator with seed: " << seed);
}

/**
 * initialize random number generator from state
 */
void ljfluid::rng(rng::gsl::gfsr4::state_type const& state)
{
    rng_.restore(state);
    LOG("restoring random number generator from state");
}

/**
 * place particles on a face-centered cubic (fcc) lattice
 */
void ljfluid::lattice()
{
    LOG("placing particles on face-centered cubic (fcc) lattice");

    // particles per 2- or 3-dimensional unit cell
    const unsigned int m = 2 * (dimension - 1);
    // lower boundary for number of particles per lattice dimension
    unsigned int n = std::pow(npart / m, 1. / dimension);
    // lower boundary for total number of lattice sites
    unsigned int N = m * std::pow(n, dimension);

    if (N < npart) {
	n += 1;
	N = m * std::pow(n, dimension);
    }
    if (N > npart) {
	LOG_WARNING("lattice not fully occupied (" << N << " sites)");
    }

    // lattice distance
    double a = box_ / n;
    // minimum distance in 2- or 3-dimensional fcc lattice
    LOG("minimum lattice distance: " << a / std::sqrt(2.));

    for (unsigned int i = 0; i < npart; ++i) {
	hvector r(a);
#ifdef DIM_3D
	// compose primitive vectors from 1-dimensional index
	r[0] *= ((i >> 2) % n) + ((i ^ (i >> 1)) & 1) / 2.;
	r[1] *= ((i >> 2) / n % n) + (i & 1) / 2.;
	r[2] *= ((i >> 2) / n / n) + (i & 2) / 4.;
#else
	// compose primitive vectors from 1-dimensional index
	r[0] *= ((i >> 1) % n) + (i & 1) / 2.;
	r[1] *= ((i >> 1) / n) + (i & 1) / 2.;
#endif
	// add particle to appropriate cell list
	compute_cell(r).push_back(particle(r, i));
	// copy position to sorted particle list
	part.r[i] = r;
    }

    // update Verlet neighbour lists
    update_neighbours();
    // reset sum over maximum velocity magnitudes to zero
    v_max_sum = 0.;
    // calculate forces, potential energy and virial equation sum
    compute_forces();
}

/**
 * set system temperature according to Maxwell-Boltzmann distribution
 */
void ljfluid::temperature(double value)
{
    LOG("initializing velocities from Maxwell-Boltzmann distribution at temperature: " << value);

    // center of mass velocity
    hvector v_cm = 0.;
    // maximum squared velocity
    double vv_max = 0.;

    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle& p, *it) {
	    // generate random Maxwell-Boltzmann distributed velocity
	    rng_.gaussian(p.v[0], p.v[1], value);
#ifdef DIM_3D
	    // Box-Muller transformation strictly generates 2 variates at once
	    rng_.gaussian(p.v[1], p.v[2], value);
#endif
	    v_cm += p.v;

	    // initialize force to zero for first leapfrog half step
	    p.f = 0.;
	}
    }

    v_cm /= npart;

    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle& p, *it) {
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
 * write parameters to HDF5 parameter group
 */
void ljfluid::attrs(H5::Group const& param) const
{
    H5xx::group node(param.createGroup("mdsim"));
    node["dimension"] = dimension;
    node["particles"] = npart;
    node["cells"] = ncell;
    node["cell_length"] = cell_length_;
    node["density"] = density_;
    node["box_length"] = box_;
    node["timestep"] = timestep_;
    node["cutoff_radius"] = r_cut;
}

/**
 * update cell lists
 */
void ljfluid::update_cells()
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
ljfluid::cell_list& ljfluid::compute_cell(hvector const& r)
{
    hvector idx = (r - floor(r / box_) * box_) / cell_length_;
    //
    // Wrapping the positional coordinates of a particle in the above way
    // has proven to reliably deliver an integer in [0, ncell - 1] after
    // conversion from floating-point for valid simulation parameters.
    //
    // However, in the extreme case of diverging potential energy (due to a
    // too-large timestep) the coordinates may jump to several orders of
    // magnitude of the box length, resulting in an out-of-bounds cell
    // index after rounding and integer conversion. Therefore, we have to
    // add integer modulo operations as a safeguard.
    //
#ifdef DIM_3D
    return cell[(unsigned int)(idx[0]) % ncell][(unsigned int)(idx[1]) % ncell][(unsigned int)(idx[2]) % ncell];
#else
    return cell[(unsigned int)(idx[0]) % ncell][(unsigned int)(idx[1]) % ncell];
#endif
}

/**
 * update neighbour lists
 */
void ljfluid::update_neighbours()
{
#ifdef DIM_3D
    for (unsigned int x = 0; x < ncell; ++x) {
	for (unsigned int y = 0; y < ncell; ++y) {
	    for (unsigned int z = 0; z < ncell; ++z) {
		foreach (particle& p, cell[x][y][z]) {
		    // empty neighbour list of particle
		    p.neighbour.clear();

		    // visit this cell
		    compute_cell_neighbours<true>(p, cell[x][y][z]);

		    // visit half of 26 neighbour cells due to pair potential
		    boost::array<cell_index, 13> neighbour = {{
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
			compute_cell_neighbours<false>(p, cell[(x + ncell + idx[0]) % ncell][(y + ncell + idx[1]) % ncell][(z + ncell + idx[2]) % ncell]);
		    }
		}
	    }
	}
    }
#else
    for (unsigned int x = 0; x < ncell; ++x) {
	for (unsigned int y = 0; y < ncell; ++y) {
	    foreach (particle& p, cell[x][y]) {
		// empty neighbour list of particle
		p.neighbour.clear();

		// visit this cell
		compute_cell_neighbours<true>(p, cell[x][y]);

		// visit half of 8 neighbour cells due to pair potential
		boost::array<cell_index, 4> neighbour = {{
		    {{  0, -1 }},
		    {{ -1, -1 }},
		    {{ -1,  0 }},
		    {{ -1, +1 }},
		}};

		// update neighbour list of particle
		foreach (cell_index const& idx, neighbour) {
		    compute_cell_neighbours<false>(p, cell[(x + ncell + idx[0]) % ncell][(y + ncell + idx[1]) % ncell]);
		}
	    }
	}
    }
#endif
}

/**
 * update neighbour list of particle
 */
template <bool same_cell>
void ljfluid::compute_cell_neighbours(particle& p, cell_list& c)
{
    foreach (particle& pp, c) {
	// skip identical particle and particle pair permutations if same cell
	if (same_cell && pp.n <= p.n)
	    continue;

	// particle distance vector
	hvector r = p.r - pp.r;
	// enforce periodic boundary conditions
	r -= round(r / box_) * box_;
	// squared particle distance
	double rr = r * r;

	// enforce cutoff radius with neighbour list skin
	if (rr >= rr_cut_skin)
	    continue;

	// add particle to neighbour list
	p.neighbour.push_back(&pp);
    }
}

/**
 * compute Lennard-Jones forces
 */
void ljfluid::compute_forces()
{
    // initialize particle forces to zero
    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle& p, *it) {
	    p.f = 0.;
	}
    }

    // potential energy
    en_pot_ = 0.;
    // virial equation sum
    virial_ = 0.;

    // iterate over all particles
    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle& p, *it) {
	    // calculate pairwise Lennard-Jones force with neighbour particles
	    foreach (particle* pp, p.neighbour) {
		// particle distance vector
		hvector r = p.r - pp->r;
		// enforce periodic boundary conditions
		r -= round(r / box_) * box_;
		// squared particle distance
		double rr = r * r;

		// enforce cutoff radius
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

    // ensure that system is still in valid state
    if (std::isinf(en_pot_)) {
	throw exception("potential energy diverged due to excessive timestep or density");
    }
}

/**
 * first leapfrog step of integration of equations of motion
 */
void ljfluid::leapfrog_half()
{
    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle& p, *it) {
	    // half step velocity
	    p.v += p.f * (timestep_ / 2.);
	    // full step position
	    p.r += p.v * timestep_;
	    // copy position to sorted particle list
	    part.r[p.n] = p.r;
	}
    }
}

/**
 * second leapfrog step of integration of equations of motion
 */
void ljfluid::leapfrog_full()
{
    // maximum squared velocity
    double vv_max = 0.;

    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	foreach (particle& p, *it) {
	    // full step velocity
	    p.v += p.f * (timestep_ / 2.);
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
void ljfluid::mdstep()
{
    // nanosecond resolution process times
    boost::array<timespec, 5> t;

    // calculate particle positions
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[0]);
    leapfrog_half();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[1]);

    if (v_max_sum * timestep_ > r_skin / 2.) {
	// update cell lists
	update_cells();
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[2]);
	// update Verlet neighbour lists
	update_neighbours();
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[3]);
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0.;

	// CPU time for update cell lists
	m_times[0] += t[2] - t[1];
	// CPU time for update Verlet neighbour lists
	m_times[1] += t[3] - t[2];
    }

    // calculate forces, potential energy and virial equation sum
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[2]);
    compute_forces();
    // calculate velocities
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[3]);
    leapfrog_full();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[4]);

    // CPU ticks for Lennard-Jones force update
    m_times[2] += t[3] - t[2];
    // CPU ticks for velocity-Verlet
    m_times[3] += (t[1] - t[0]) + (t[4] - t[3]);
    // CPU ticks for MD simulation step
    m_times[4] += t[4] - t[0];
}

/**
 * sample trajectory
 */
void ljfluid::sample(mdsim_sample_visitor visitor) const
{
    visitor(part.r, part.v, en_pot_, virial_);
}

/*
 * returns and resets CPU tick statistics
 */
perf_counters ljfluid::times()
{
    perf_counters times(m_times);
    // reset performance counters
    for (unsigned int i = 0; i < m_times.size(); ++i) {
	m_times[i].clear();
    }
    return times;
}

} // namespace mdsim
