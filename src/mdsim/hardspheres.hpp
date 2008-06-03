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
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>
#include "H5param.hpp"
#include "exception.hpp"
#include "gsl_rng.hpp"
#include "log.hpp"


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
    /**
     * particle state
     */
    struct particle
    {
	/** particle position immediately after most recently processed event for this particle */
	T r;
	/** particle velocity immediately after most recently processed event for this particle */
	T v;
	/** time of that event */
	double t;
	/** event counter */
	uint64_t count;

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
	    /** reexamine due to periodic boundary conditions */
	    BOUNDARY,
	} type;

	/** collision event partner */
	unsigned int n2;
	/** copy of event counter of partner at time of event */
	uint64_t count2;
    };

    /** particle event queue item with event time and particle */
    typedef std::pair<double, unsigned int> event_queue_item;

public:
    /** initialize current simulation time */
    hardspheres() : time_(0.) {}

    /** set number of particles */
    void particles(unsigned int value);
    /** set pair separation at which particle collision occurs */
    void pair_separation(double value);
    /** set particle density */
    void density(double value);
    /** set periodic box length */
    void box(double value);

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

    /** advance phase space state to given sample time */
    void mdstep(double const& time);
    /** sample trajectory */
    template <typename V> void sample(V visitor) const;

private:
    /** schedule next particle event starting at given time */
    void schedule_event(unsigned int const& n, double const& t);
    /** process next particle event */
    void process_event(unsigned int const& n, double const& t);

private:
    /** number of particles */
    unsigned int npart;
    /** pair separation at which particle collision occurs */
    double pair_sep_;
    /** particle density */
    double density_;
    /** periodic box length */
    double box_;

    /** particle states */
    std::vector<particle> part;
    /** particle event list with next event for each particle */
    std::vector<event> event_list;
    /** time-ordered particle event queue */
    std::priority_queue<event_queue_item, std::vector<event_queue_item>, std::greater<event_queue_item> > event_queue;
    /** current simulation time */
    double time_;

    /** particle positions at sample time */
    std::vector<T> r_;
    /** particle velocities at sample time */
    std::vector<T> v_;

    /** random number generator */
    rng::gsl::gfsr4 rng_;
    /** squared pair separation */
    double pair_sep_sq;
    /** periodic boundary condition reexamination event distance */
    double bound_dist;
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

    // derive periodic boundary condition reexamination event distance
    bound_dist = (box_ / 2. - pair_sep_) / 2.;
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

    // derive periodic boundary condition reexamination event distance
    bound_dist = (box_ / 2. - pair_sep_) / 2.;
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

    for (unsigned int i = 0; i < npart; ++i) {
	// set particle position at simulation time zero
	part[i].r = r_[i];
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
	// set particle time
	part[i].t = 0.;
	// copy particle position at sample time zero
	r_[i] = part[i].r;
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
	schedule_event(i, 0.);
    }
}

/**
 * schedule next particle event starting at given time
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::schedule_event(unsigned int const& n, double const& t)
{
    // time interval till collision
    double dt = std::numeric_limits<double>::max();
    // collision partner
    int n2 = -1;

    // find earliest particle collision event starting at given time
    for (unsigned int j = 0; j < npart; ++j) {
	// skip same particle
	if (j == n)
	    continue;

	// particle distance vector at given time
	T dr = part[j].r + part[j].v * (t - part[j].t) - (part[n].r + part[n].v * (t - part[n].t));
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

    // check if collision event has been found
    if (n2 >= 0) {
	// generate particle collision event
	event_list[n].type = event::COLLISION;
	event_list[n].t = t + dt;
	event_list[n].n2 = n2;
	event_list[n].count2 = part[n2].count;
    }
    else {
	// generate periodic boundary condition reexamination event
	event_list[n].type = event::BOUNDARY;
	// minimum time required to traverse boundary event distance in any direction
	dt = std::min(bound_dist / fabs(part[n].v.x), bound_dist / fabs(part[n].v.y));
#ifdef DIM_3D
	dt = std::min(dt, bound_dist / fabs(part[n].v.z));
#endif
	event_list[n].t = t + dt;
    }

    // schedule particle event
    event_queue.push(event_queue_item(event_list[n].t, n));
}

/*
 * process next particle event
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::process_event(unsigned int const& n, double const& t)
{
    if (event_list[n].t != t)
	// discard invalidated event
	return;

    if (event_list[n].type == event::COLLISION) {
	// collision partner particle number
	const unsigned int n2 = event_list[n].n2;

	// check if partner participated in another collision before this event
	if (part[n2].count != event_list[n].count2) {
	    // schedule next event for this particle
	    schedule_event(n, t);
	    return;
	}

	// update positions to current simulation time
	part[n].r += part[n].v * (t - part[n].t);
	part[n2].r += part[n2].v * (t - part[n2].t);
	// update particle times
	part[n].t = t;
	part[n2].t = t;

	// particle distance vector
	T dr = part[n2].r - part[n].r;
	// enforce periodic boundary conditions
	dr -= round(dr / box_) * box_;
	// velocity difference before collision
	T dv = part[n].v - part[n2].v;
	// velocity difference after collision without dissipation
	dv = dr * (dr * dv) / (dr * dr);

	// update velocities to current simulation time
	part[n].v -= dv;
	part[n2].v += dv;
	// update particle event counters
	part[n].count++;
	part[n2].count++;

	// schedule next event for each particle
	schedule_event(n, t);
	schedule_event(n2, t);
    }
    else if (event_list[n].type == event::BOUNDARY) {
	// update particle position to current simulation time
	part[n].r += part[n].v * (t - part[n].t);
	// update particle time
	part[n].t = t;
	// update particle event counter
	part[n].count++;

	// schedule next event for particle
	schedule_event(n, t);
    }
}

/**
 * advance phase space state to given sample time
 */
template <unsigned dimension, typename T>
void hardspheres<dimension, T>::mdstep(double const& sample_time)
{
    // process particle event queue till sample time
    while (event_queue.top().first <= sample_time) {
	if (event_queue.top().first < time_) {
	    throw exception("simulation time is running backwards");
	}

	// update current simulation time to particle event time
	time_ = event_queue.top().first;

	process_event(event_queue.top().second, time_);
	event_queue.pop();
    }

    // sample phase space at given time
    for (unsigned int i = 0; i < npart; ++i) {
	// particle position
	r_[i] = part[i].r + part[i].v * (sample_time - part[i].t);
	// enforce periodic boundary conditions
	r_[i] -= floor(r_[i] / box_) * box_;
	// particle velocity
	v_[i] = part[i].v;
    }
}

/**
 * sample trajectory
 */
template <unsigned dimension, typename T>
template <typename V>
void hardspheres<dimension, T>::sample(V visitor) const
{
    visitor(r_, v_);
}

} // namespace mdsim

#undef foreach

#endif /* ! MDSIM_HARDSPHERES_HPP */
