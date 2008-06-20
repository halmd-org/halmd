/* Hard Spheres simulation
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

#ifndef MDSIM_HARDSPHERES_HPP
#define MDSIM_HARDSPHERES_HPP

#include <algorithm>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <queue>
#include <vector>
#include "H5param.hpp"
#include "exception.hpp"
#include "gsl_rng.hpp"
#include "log.hpp"
#include "perf.hpp"
#include "time.hpp"


#define foreach BOOST_FOREACH

namespace mdsim
{

/**
 * Hard Spheres simulation
 */
template <unsigned dimension, typename T>
class hardspheres
{
    //
    // Details of the implementation are described in
    //
    // S. Miller, S. Luding,
    // Event-driven molecular dynamics in parallel,
    // Journal of Computational Physics 193 (2003) 306-316
    //
    // M. P. Allen, D. Frenkel & J. Talbot,
    // Molecular dynamics simulation using hard particles,
    // Computer Physics Reports 9 (1989) 301-353
    //

public:
    typedef std::list<unsigned int> cell_type;
    typedef boost::array<unsigned int, dimension> cell_index;

    /**
     * particle state
     */
    struct particle
    {
	/** periodically reduced particle position */
	T r;
	/** periodically extended particle position */
	T R;
	/** particle velocity */
	T v;
	/** time of that event */
	double t;
	/** event counter */
	uint64_t count;
	/** cell which particle belongs to */
	cell_index cell;

	/** initialize event counter to zero */
	particle() : count(0) {}
    };

    /**
     * particle event list item
     */
    struct event
    {
	/** time of event */
	double t;
	/** event type */
	enum {
	    /** collision with other particle */
	    COLLISION,
	    /** cell boundary */
	    CELL,
	} type;

	/** collision event partner */
	unsigned int n2;
	/** cell boundary */
	cell_index cell2;
	/** copy of event counter of partner at time of event */
	uint64_t count2;
    };

    /** particle event queue item with event time and particle */
    typedef std::pair<double, unsigned int> event_queue_item;

public:
    /** set number of particles */
    void particles(unsigned int value);
    /** set pair separation at which particle collision occurs */
    void pair_separation(double value);
    /** set particle density */
    void density(double value);
    /** set periodic box length */
    void box(double value);
    /** initialize cells */
    void init_cell();

    /** set system state from phase space sample */
    template <typename V> void restore(V visitor);
    /** initialize random number generator with seed */
    void rng(unsigned int seed);
    /** initialize random number generator from state */
    void rng(rng::gsl::gfsr4::state_type const& state);
    /** place particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(double value);
    /** initialize event list */
    void init_event_list();

    /** returns number of particles */
    unsigned int const& particles() const;
    /** returns pair separation at which particle collision occurs */
    double const& pair_separation() const;
    /** returns particle density */
    double const& density() const;
    /** returns periodic box length */
    double const& box() const;
    /** returns number of cells per dimension */
    unsigned int const& cells() const;
    /** returns cell length */
    double const& cell_length();

    /** advance phase space state to given sample time */
    void mdstep(double sample_time);
    /** sample trajectory */
    template <typename V> void sample(V visitor) const;
    /** get execution time statistics */
    perf_type const& times() const;

private:
    /** schedule next particle event starting at given time */
    void schedule_event(unsigned int n);
    /** process particle collision event */
    void process_collision_event(unsigned int n);
    /** process cell boundary event */
    void process_cell_event(unsigned int n);
    /** returns cell which a particle belongs to */
    cell_index compute_cell(T const& r);
    /** compute next collision event with particles of given cell starting at given time within given time interval */
    void compute_collision_event(unsigned int n, cell_type const& cell);
    /** compute next cell boundary event starting at given time within given time interval */
    void compute_cell_event(unsigned int n);

private:
    /** number of particles */
    unsigned int npart;
    /** pair separation at which particle collision occurs */
    double pair_sep_;
    /** particle density */
    double density_;
    /** periodic box length */
    double box_;
    /** number of cells per dimension */
    unsigned int ncell;
    /** cell length */
    double cell_length_;

    /** particle states */
    std::vector<particle> part;
    /** cells */
    boost::multi_array<cell_type, dimension> cell_;
    /** particle event list with next event for each particle */
    std::vector<event> event_list;
    /** time-ordered particle event queue */
    std::priority_queue<event_queue_item, std::vector<event_queue_item>, std::greater<event_queue_item> > event_queue;

    /** periodically reduced particle positions at sample time */
    std::vector<T> r_;
    /** periodically extended particle positions at sample time */
    std::vector<T> R_;
    /** particle velocities at sample time */
    std::vector<T> v_;
    /** impulsive limit of the virial expression sum */
    double virial_;

    /** random number generator */
    rng::gsl::gfsr4 rng_;
    /** squared pair separation */
    double pair_sep_sq;

    /** execution time statistics */
    perf_type times_;
};

/**
 * set number of particles in system
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::particles(unsigned int value)
{
    if (value < 1) {
	throw exception("number of particles must be non-zero");
    }
    npart = value;
    LOG("number of particles: " << npart);

    try {
	part.resize(npart);
	r_.resize(npart);
	R_.resize(npart);
	v_.resize(npart);
    }
    catch (std::bad_alloc const&) {
	throw exception("failed to allocate particle states");
    }
}

/**
 * set pair separation at which particle collision occurs
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::pair_separation(double value)
{
    if (value <= 0.) {
	throw exception("pair separation must be greater than zero");
    }
    pair_sep_ = value;
    LOG("pair separation: " << pair_sep_);

    // squared pair separation
    pair_sep_sq = pair_sep_ * pair_sep_;
}

/**
 * set particle density
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::density(double value)
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
void hardspheres<dimension, T>::box(double value)
{
    box_ = value;
    LOG("periodic box length: " << box_);

    // derive particle density
    density_ = npart / std::pow(box_, 1. * dimension);
    LOG("particle density: " << density_);
}

/**
 * initialize cells
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::init_cell()
{
#ifdef DIM_3D
    // FIXME optimal number of cells
    ncell = std::min(cbrt(npart * 8.), std::floor(box_ / pair_sep_));
#else
    ncell = std::min(sqrt(npart * 1.5), std::floor(box_ / pair_sep_));
#endif
    LOG("number of cells per dimension: " << ncell);

    if (ncell < 3) {
	throw exception("number of cells per dimension must be at least 3");
    }

    try {
#ifdef DIM_3D
	cell_.resize(boost::extents[ncell][ncell][ncell]);
#else
	cell_.resize(boost::extents[ncell][ncell]);
#endif
    }
    catch (std::bad_alloc const&) {
	throw exception("failed to allocate cells");
    }

    // derive cell length
    cell_length_ = box_ / ncell;
    LOG("cell length: " << cell_length_);
}

/**
 * set system state from phase space sample
 */
template <unsigned dimension, typename T>
template <typename V>
void hardspheres<dimension, T>::restore(V visitor)
{
    // copy particle positions and velocities at sample time zero
    visitor(r_, v_);
    // replicate to periodically extended coordinates
    std::copy(r_.begin(), r_.end(), R_.begin());

    for (unsigned int i = 0; i < npart; ++i) {
	// set periodically reduced particle position at simulation time zero
	part[i].r = r_[i];
	// set periodically extended particle position at simulation time zero
	part[i].R = R_[i];
	// set cell which particle belongs to
	part[i].cell = compute_cell(part[i].r);
	// add particle to cell
	cell_(part[i].cell).push_back(i);
	// set particle velocity at simulation time zero
	part[i].v = v_[i];
	// set particle time
	part[i].t = 0.;
    }
}

/**
 * initialize random number generator with seed
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::rng(unsigned int seed)
{
    rng_.set(seed);
    LOG("initializing random number generator with seed: " << seed);
}

/**
 * initialize random number generator from state
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::rng(rng::gsl::gfsr4::state_type const& state)
{
    rng_.restore(state);
    LOG("restoring random number generator from state");
}

/**
 * place particles on a face-centered cubic (fcc) lattice
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::lattice()
{
    LOG("placing particles on face-centered cubic (fcc) lattice");

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
#ifdef DIM_3D
	// compose primitive vectors from 1-dimensional index
	part[i].r.x = ((i >> 2) % n) + ((i ^ (i >> 1)) & 1) / 2.;
	part[i].r.y = ((i >> 2) / n % n) + (i & 1) / 2.;
	part[i].r.z = ((i >> 2) / n / n) + (i & 2) / 4.;
#else
	// compart[i]ose primitive vectors from 1-dimensional index
	part[i].r.x = ((i >> 1) % n) + (i & 1) / 2.;
	part[i].r.y = ((i >> 1) / n) + (i & 1) / 2.;
#endif
	// scale by lattice distance
	part[i].r *= a;
	// set periodically extended particle position
	part[i].R = part[i].r;
	// set cell which particle belongs to
	part[i].cell = compute_cell(part[i].r);
	// add particle to cell
	cell_(part[i].cell).push_back(i);
	// set particle time
	part[i].t = 0.;
	// copy periodically reduced particle position at sample time zero
	r_[i] = part[i].r;
	// copy periodically extended particle position at sample time zero
	R_[i] = part[i].R;
    }
}

/**
 * set system temperature according to Maxwell-Boltzmann distribution
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::temperature(double value)
{
    LOG("initializing velocities from Maxwell-Boltzmann distribution at temperature: " << value);

    // center of mass velocity
    T v_cm = 0.;

    foreach (particle& p, part) {
	// generate random Maxwell-Boltzmann distributed velocity
	rng_.gaussian(p.v.x, p.v.y, value);
#ifdef DIM_3D
	// Box-Muller transformation strictly generates 2 variates at once
	rng_.gaussian(p.v.y, p.v.z, value);
#endif
	v_cm += p.v;
    }

    v_cm /= npart;

    for (unsigned int i = 0; i < npart; ++i) {
	// set center of mass velocity to zero
	part[i].v -= v_cm;
	// copy particle velocity at sample time zero
	v_[i] = part[i].v;
    }
}

/**
 * returns number of particles
 */
template <unsigned dimension, typename T>
unsigned int const& hardspheres<dimension, T>::particles() const
{
    return npart;
}

/**
 * returns pair separation at which particle collision occurs
 */
template <unsigned dimension, typename T>
double const& hardspheres<dimension, T>::pair_separation() const
{
    return pair_sep_;
}

/**
 * returns particle density
 */
template <unsigned dimension, typename T>
double const& hardspheres<dimension, T>::density() const
{
    return density_;
}

/**
 * returns periodic box length
 */
template <unsigned dimension, typename T>
double const& hardspheres<dimension, T>::box() const
{
    return box_;
}

/**
 * returns number of cells per dimension
 */
template <unsigned dimension, typename T>
unsigned int const& hardspheres<dimension, T>::cells() const
{
    return ncell;
}

/**
 * returns cell length
 */
template <unsigned dimension, typename T>
double const& hardspheres<dimension, T>::cell_length()
{
    return cell_length_;
}

/**
 * initialize event list
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::init_event_list()
{
    try {
	event_list.resize(npart);
    }
    catch (std::bad_alloc const&) {
	throw exception("failed to allocate event list");
    }

    // schedule next event for each particle
    for (unsigned int i = 0; i < npart; ++i) {
	schedule_event(i);
    }
}

/**
 * compute next collision event with particles of given cell
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::compute_collision_event(const unsigned int n, cell_type const& cell)
{
    double dt = std::numeric_limits<double>::max();
    int n2 = -1;

    // iterate over particles in cell
    foreach (unsigned int j, cell) {
	// skip same particle if in same cell
	if (j == n)
	    continue;

	// particle distance vector at time of first particle
	T dr = part[j].r + part[j].v * (part[n].t - part[j].t) - part[n].r;
	// enforce periodic boundary conditions
	dr -= round(dr / box_) * box_;
	// velocity difference at given time
	T dv = part[j].v - part[n].v;

	// check particle collision constraint
	const double drdv = dr * dv;
	if (drdv >= 0.)
	    // no particle collision in future
	    continue;
	const double dvdv = dv * dv;
	const double rad = (drdv * drdv) - dvdv * ((dr * dr) - pair_sep_sq);
	if (rad < 0.)
	    // no particle collision in future
	    continue;
	const double dt_ = (- drdv - std::sqrt(rad)) / dvdv;
	if (dt_ < 0.)
	    // no particle collision in future
	    continue;

	// particles will collide in the future in reference to given time
	if (dt_ < dt) {
	    // set smallest collision time interval
	    dt = dt_;
	    // set partner participating in that collision
	    n2 = j;
	}
    }

    if (n2 < 0)
	// no collision with particles in cell
	return;

    if (dt < event_list[n].t - part[n].t) {
	// generate particle collision event
	event_list[n].type = event::COLLISION;
	event_list[n].t = part[n].t + dt;
	event_list[n].n2 = n2;
	event_list[n].count2 = part[n2].count;
    }
}

/**
 * compute next cell boundary event
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::compute_cell_event(const unsigned int n)
{
    T dt3(std::numeric_limits<double>::max());
    cell_index cell2 = part[n].cell;

    if (part[n].v.x < 0.) {
	dt3.x = (part[n].cell[0] * cell_length_ - part[n].r.x) / part[n].v.x;
	cell2[0] = (cell2[0] + ncell - 1) % ncell;
    }
    else if (part[n].v.x > 0.) {
	dt3.x = ((part[n].cell[0] + 1) * cell_length_ - part[n].r.x) / part[n].v.x;
	cell2[0] = (cell2[0] + 1) % ncell;
    }
    if (part[n].v.y < 0.) {
	dt3.y = (part[n].cell[1] * cell_length_ - part[n].r.y) / part[n].v.y;
	cell2[1] = (cell2[1] + ncell - 1) % ncell;
    }
    else if (part[n].v.y > 0.) {
	dt3.y = ((part[n].cell[1] + 1) * cell_length_ - part[n].r.y) / part[n].v.y;
	cell2[1] = (cell2[1] + 1) % ncell;
    }
#ifdef DIM_3D
    if (part[n].v.z < 0.) {
	dt3.z = (part[n].cell[2] * cell_length_ - part[n].r.z) / part[n].v.z;
	cell2[2] = (cell2[2] + ncell - 1) % ncell;
    }
    else if (part[n].v.z > 0.) {
	dt3.z = ((part[n].cell[2] + 1) * cell_length_ - part[n].r.z) / part[n].v.z;
	cell2[2] = (cell2[2] + 1) % ncell;
    }
#endif

#ifdef DIM_3D
    const double dt = std::min(std::min(dt3.x, dt3.y), dt3.z);
#else
    const double dt = std::min(dt3.x, dt3.y);
#endif

    if (dt < event_list[n].t - part[n].t) {
	// generate cell boundary event
	event_list[n].t = part[n].t + dt;
	event_list[n].type = event::CELL;
	event_list[n].cell2[0] = (dt3.x == dt) ? cell2[0] : part[n].cell[0];
	event_list[n].cell2[1] = (dt3.y == dt) ? cell2[1] : part[n].cell[1];
#ifdef DIM_3D
	event_list[n].cell2[2] = (dt3.z == dt) ? cell2[2] : part[n].cell[2];
#endif
    }
}

/**
 * schedule next particle event starting at given time
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::schedule_event(const unsigned int n)
{
    // upper boundary for time of next particle event
    event_list[n].t = std::numeric_limits<double>::max();

    // compute next cell boundary event
    compute_cell_event(n);
    // compute next collision event with particles of this cell
    compute_collision_event(n, cell_(part[n].cell));

#ifdef DIM_3D
    // compute next collision event with particles of neighbour cells
    const boost::array<cell_index, 26> neighbour = {{
	{{  0, -1,  0 }},
	{{  0, +1,  0 }},
	{{ -1, -1,  0 }},
	{{ -1,  0,  0 }},
	{{ -1, +1,  0 }},
	{{ +1, -1,  0 }},
	{{ +1,  0,  0 }},
	{{ +1, +1,  0 }},
	{{  0, -1, -1 }},
	{{  0, +1, -1 }},
	{{  0, +1, +1 }},
	{{ -1, -1, -1 }},
	{{ -1,  0, -1 }},
	{{ -1, +1, -1 }},
	{{ +1, -1, -1 }},
	{{ +1,  0, -1 }},
	{{ +1, +1, -1 }},
	{{  0, -1, +1 }},
	{{ -1, -1, +1 }},
	{{ -1,  0, +1 }},
	{{ -1, +1, +1 }},
	{{ +1, -1, +1 }},
	{{ +1,  0, +1 }},
	{{ +1, +1, +1 }},
	{{  0,  0, -1 }},
	{{  0,  0, +1 }},
    }};

    foreach (cell_index const& idx, neighbour) {
	compute_collision_event(n, cell_[(part[n].cell[0] + ncell + idx[0]) % ncell][(part[n].cell[1] + ncell + idx[1]) % ncell][(part[n].cell[2] + ncell + idx[2]) % ncell]);
    }
#else
    // compute next collision event with particles of neighbour cells
    const boost::array<cell_index, 8> neighbour = {{
	{{  0, -1 }},
	{{ -1, -1 }},
	{{ -1,  0 }},
	{{ -1, +1 }},
	{{  0, +1 }},
	{{ +1, -1 }},
	{{ +1,  0 }},
	{{ +1, +1 }},
    }};

    foreach (cell_index const& idx, neighbour) {
	compute_collision_event(n, cell_[(part[n].cell[0] + ncell + idx[0]) % ncell][(part[n].cell[1] + ncell + idx[1]) % ncell]);
    }
#endif

    // schedule particle event
    event_queue.push(event_queue_item(event_list[n].t, n));
}

/*
 * process particle collision event
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::process_collision_event(const unsigned int n1)
{
    const T dr1 = part[n1].v * (event_list[n1].t - part[n1].t);
    // update periodically extended particle position
    part[n1].R += dr1;
    // update periodically reduced particle position to given time
    part[n1].r += dr1;
    // update particle time
    part[n1].t = event_list[n1].t;

    // collision partner particle number
    const unsigned int n2 = event_list[n1].n2;

    // check if partner participated in another collision before this event
    if (part[n2].count != event_list[n1].count2) {
	// schedule next event for this particle
	schedule_event(n1);
	return;
    }

    const T dr2 = part[n2].v * (event_list[n1].t - part[n2].t);
    // update periodically extended particle position
    part[n2].R += dr2;
    // update periodically reduced particle position to given time
    part[n2].r += dr2;
    // update particle time
    part[n2].t = event_list[n1].t;

    // particle distance vector
    T dr = part[n2].r - part[n1].r;
    // enforce periodic boundary conditions
    dr -= round(dr / box_) * box_;
    // velocity difference before collision
    T dv = part[n1].v - part[n2].v;
    // velocity difference after collision without dissipation
    dv = dr * (dr * dv) / (dr * dr);

    // update velocities to current simulation time
    part[n1].v -= dv;
    part[n2].v += dv;

    // add contribution to impulsive limit of the virial expression sum
    virial_ += dr * dv;

    // update particle event counters
    part[n1].count++;
    part[n2].count++;

    // schedule next event for each particle
    schedule_event(n1);
    schedule_event(n2);
}

/*
 * process cell boundary event
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::process_cell_event(const unsigned int n)
{
    const T dr = part[n].v * (event_list[n].t - part[n].t);
    // update periodically extended particle position
    part[n].R += dr;
    // update periodically reduced particle position to given time
    part[n].r += dr;
    // enforce periodic boundary conditions
    if (part[n].cell[0] == ncell - 1 && event_list[n].cell2[0] == 0)
	part[n].r.x -= box_;
    if (part[n].cell[0] == 0 && event_list[n].cell2[0] == ncell - 1)
	part[n].r.x += box_;
    if (part[n].cell[1] == ncell - 1 && event_list[n].cell2[1] == 0)
	part[n].r.y -= box_;
    if (part[n].cell[1] == 0 && event_list[n].cell2[1] == ncell - 1)
	part[n].r.y += box_;
#ifdef DIM_3D
    if (part[n].cell[2] == ncell - 1 && event_list[n].cell2[2] == 0)
	part[n].r.z -= box_;
    if (part[n].cell[2] == 0 && event_list[n].cell2[2] == ncell - 1)
	part[n].r.z += box_;
#endif
    // update particle time
    part[n].t = event_list[n].t;

    // remove particle from old cell
    cell_(part[n].cell).remove(n);
    // update particle cell
    part[n].cell = event_list[n].cell2;
    // add particle to cell
    cell_(part[n].cell).push_back(n);

    // schedule next event for particle
    schedule_event(n);
}

/**
 * returns cell which a particle belongs to
 */
template <unsigned dimension, typename T>
typename hardspheres<dimension, T>::cell_index hardspheres<dimension, T>::compute_cell(T const& r)
{
    T cellf = r / cell_length_;
#ifdef DIM_3D
    cell_index cell = {{ int(cellf.x), int(cellf.y), int(cellf.z) }};
#else
    cell_index cell = {{ int(cellf.x), int(cellf.y) }};
#endif
    return cell;
}

/**
 * advance phase space state to given sample time
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::mdstep(const double sample_time)
{
    boost::array<timer, 3> t;
    t[0].start();

    // impulsive limit of the virial expression sum
    virial_ = 0.;

    // process particle event queue till sample time
    t[1].start();
    while (event_queue.top().first <= sample_time) {
	if (event_queue.top().first != event_list[event_queue.top().second].t) {
	    // discard invalidated event
	    event_queue.pop();
	    continue;
	}

	switch (event_list[event_queue.top().second].type) {
	  case event::COLLISION:
	    // process particle collision event
	    process_collision_event(event_queue.top().second);
	    break;

	  case event::CELL:
	    // process cell boundary event
	    process_cell_event(event_queue.top().second);
	    break;
	}
	event_queue.pop();
    }
    t[1].stop();
    times_["host"]["queue"] += t[1].elapsed();

    virial_ /= npart;

    // sample phase space at given time
    t[2].start();
    for (unsigned int i = 0; i < npart; ++i) {
	const T dr = part[i].v * (sample_time - part[i].t);
	// periodically extended particle position
	R_[i] = part[i].R + dr;
	// periodically reduced particle position
	r_[i] = part[i].r + dr;
	// enforce periodic boundary conditions
	r_[i] -= floor(r_[i] / box_) * box_;
	// particle velocity
	v_[i] = part[i].v;
    }
    t[2].stop();
    times_["host"]["sample"] += t[2].elapsed();

    t[0].stop();
    times_["host"]["mdstep"] += t[0].elapsed();
}

/**
 * sample trajectory
 */
template <unsigned dimension, typename T>
template <typename V>
void hardspheres<dimension, T>::sample(V visitor) const
{
    visitor(r_, R_, v_, virial_);
}
 
/*
 * get execution time statistics
 */
template <unsigned dimension, typename T>
perf_type const& hardspheres<dimension, T>::times() const
{
    return times_;
}

} // namespace mdsim

#undef foreach

#endif /* ! MDSIM_HARDSPHERES_HPP */
