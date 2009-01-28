/* Hard Spheres simulation
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_MDSIM_HARDSPHERE_HPP
#define LJGPU_MDSIM_HARDSPHERE_HPP

#include <algorithm>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <ljgpu/mdsim/base.hpp>
#include <ljgpu/mdsim/traits.hpp>
#include <ljgpu/rng/gsl_rng.hpp>
#include <ljgpu/sample/perf.hpp>
#include <ljgpu/util/exception.hpp>
#include <ljgpu/util/log.hpp>
#include <ljgpu/util/timer.hpp>
#include <queue>
#include <sys/times.h>
#include <vector>

namespace ljgpu
{

template <typename mdsim_impl>
class hardsphere;

/**
 * Hard Spheres simulation
 */
template<int dimension>
class hardsphere<hardsphere_impl<dimension> > : public mdsim_base<hardsphere_impl<dimension> >
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
    typedef mdsim_base<hardsphere_impl<dimension> > _Base;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::sample_type sample_type;
    typedef typename _Base::host_sample_type host_sample_type;
    typedef typename _Base::energy_sample_type energy_sample_type;
    typedef typename sample_type::sample_visitor sample_visitor;

    typedef std::list<unsigned int> cell_type;
    typedef boost::array<unsigned int, dimension> cell_index;

    /**
     * particle state
     */
    struct particle
    {
	/** periodically reduced particle position */
	vector_type r;
	/** periodically extended particle position */
	vector_type R;
	/** particle velocity */
	vector_type v;
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
    template <typename T>
    void particles(T const& value);
    /** set pair separation at which particle collision occurs */
    void pair_separation(double value);
    /** initialize cells */
    void init_cells();
    /** set simulation timestep */
    void timestep(double value);

    /** set system state from phase space sample */
    void state(sample_visitor visitor);
    /** initialize random number generator with seed */
    void rng(unsigned int seed);
    /** initialize random number generator from state */
    void rng(gsl::gfsr4::state_type const& state);
    /** place particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(double value);

    /** returns number of particles */
    unsigned int const& particles() const { return npart; }
    /** returns trajectory sample */
    sample_type const& sample() const { return m_sample; }
    /** returns pair separation at which particle collision occurs */
    double const& pair_separation() const { return pair_sep_; }
    /** returns number of cells per dimension */
    unsigned int const& cells() const { return ncell; }
    /** returns cell length */
    double const& cell_length() { return cell_length_; }
    /** returns simulation timestep */
    double const& timestep() { return timestep_; }

    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

    /** initialize event list */
    void init_event_list();
    /** advance phase space state to given sample time */
    void mdstep();
    /** sample phase space on host */
    void sample(host_sample_type& sample) const;
    /** sample thermodynamic equilibrium properties */
    void sample(energy_sample_type& sample) const;

private:
    /** schedule next particle event starting at given time */
    void schedule_event(unsigned int n);
    /** process particle collision event */
    void process_collision_event(unsigned int n);
    /** process cell boundary event */
    void process_cell_event(unsigned int n);
    /** returns cell which a particle belongs to */
    cell_index compute_cell(vector_type const& r);
    /** compute next collision event with particles of given cell starting at given time within given time interval */
    void compute_collision_event(unsigned int n, cell_type const& cell);
    /** compute next cell boundary event starting at given time within given time interval */
    void compute_cell_event(unsigned int n);

private:
    using _Base::npart;
    using _Base::box_;
    using _Base::density_;
    using _Base::m_sample;
    using _Base::m_times;

    /** pair separation at which particle collision occurs */
    double pair_sep_;
    /** number of cells per dimension */
    unsigned int ncell;
    /** cell length */
    double cell_length_;
    /** simulation timestep */
    double timestep_;

    /** particle states */
    std::vector<particle> part;
    /** cells */
    boost::multi_array<cell_type, dimension> cell_;
    /** particle event list with next event for each particle */
    std::vector<event> event_list;
    /** time-ordered particle event queue */
    std::priority_queue<event_queue_item, std::vector<event_queue_item>, std::greater<event_queue_item> > event_queue;
    /** current simulation step */
    uint64_t step_;
    /** virial equation sum per particle */
    double virial;

    /** random number generator */
    gsl::gfsr4 rng_;
    /** squared pair separation */
    double pair_sep_sq;
};

/**
 * set number of particles in system
 */
template <int dimension>
template <typename T>
void hardsphere<hardsphere_impl<dimension> >::particles(T const& value)
{
    _Base::particles(value);

    try {
	part.resize(npart);
	m_sample[0].r.resize(npart);
	m_sample[0].v.resize(npart);
    }
    catch (std::bad_alloc const&) {
	throw exception("failed to allocate particle states");
    }
}

/**
 * set pair separation at which particle collision occurs
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::pair_separation(double value)
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
 * initialize cells
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::init_cells()
{
    // FIXME optimal number of cells
    if (dimension == 3)
	ncell = std::min(cbrt(npart * 8.), std::floor(box_ / pair_sep_));
    else
	ncell = std::min(sqrt(npart * 1.5), std::floor(box_ / pair_sep_));
    LOG("number of cells per dimension: " << ncell);

    if (ncell < 3) {
	throw exception("number of cells per dimension must be at least 3");
    }

    try {
	cell_index size;
	std::fill(size.begin(), size.end(), ncell);
	cell_.resize(size);
    }
    catch (std::bad_alloc const&) {
	throw exception("failed to allocate cells");
    }

    // derive cell length
    cell_length_ = box_ / ncell;
    LOG("cell length: " << cell_length_);
}

/**
 * set simulation timestep
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::timestep(double value)
{
    timestep_ = value;
    LOG("simulation timestep: " << timestep_);
}

/**
 * set system state from phase space sample
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::state(sample_visitor visitor)
{
    _Base::state(visitor);

    for (unsigned int i = 0; i < npart; ++i) {
	// set periodically reduced particle position at simulation time zero
	part[i].r = m_sample[0].r[i];
	// set periodically extended particle position at simulation time zero
	part[i].R = part[i].r;
	// set cell which particle belongs to
	part[i].cell = compute_cell(part[i].r);
	// add particle to cell
	cell_(part[i].cell).push_back(i);
	// set particle velocity at simulation time zero
	part[i].v = m_sample[0].v[i];
	// set particle time
	part[i].t = 0.;
    }
}

/**
 * initialize random number generator with seed
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::rng(unsigned int seed)
{
    rng_.set(seed);
    LOG("initializing random number generator with seed: " << seed);
}

/**
 * initialize random number generator from state
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::rng(gsl::gfsr4::state_type const& state)
{
    rng_.restore(state);
    LOG("restoring random number generator from state");
}

/**
 * place particles on a face-centered cubic (fcc) lattice
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::lattice()
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
    const double a = box_ / n;
    // minimum distance in 2- or 3-dimensional fcc lattice
    const double dist = a / std::sqrt(2.);
    LOG("minimum lattice distance: " << dist);

    // ensure that particles do not overlap
    if (dist < pair_sep_) {
	throw exception("minimum lattice distance smaller than pair separation");
    }

    for (unsigned int i = 0; i < npart; ++i) {
	// compose primitive vectors from 1-dimensional index
	if (dimension == 3) {
	    part[i].r[0] = ((i >> 2) % n) + ((i ^ (i >> 1)) & 1) / 2.;
	    part[i].r[1] = ((i >> 2) / n % n) + (i & 1) / 2.;
	    part[i].r[2] = ((i >> 2) / n / n) + (i & 2) / 4.;
	}
	else {
	    part[i].r[0] = ((i >> 1) % n) + (i & 1) / 2.;
	    part[i].r[1] = ((i >> 1) / n) + (i & 1) / 2.;
	}
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
    }
}

/**
 * set system temperature according to Maxwell-Boltzmann distribution
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::temperature(double value)
{
    LOG("initializing velocities from Maxwell-Boltzmann distribution at temperature: " << value);

    // center of mass velocity
    vector_type v_cm = 0.;

    BOOST_FOREACH(particle& p, part) {
	// generate random Maxwell-Boltzmann distributed velocity
	rng_.gaussian(p.v[0], p.v[1], value);
	if (dimension == 3) {
	    // Box-Muller transformation strictly generates 2 variates at once
	    rng_.gaussian(p.v[1], p.v[2], value);
	}
	v_cm += p.v;
    }

    v_cm /= npart;

    for (unsigned int i = 0; i < npart; ++i) {
	// set center of mass velocity to zero
	part[i].v -= v_cm;
    }
}

/**
 * write parameters to HDF5 parameter group
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::param(H5param& param) const
{
    _Base::param(param);

    H5xx::group node(param["mdsim"]);
    node["pair_separation"] = pair_sep_;
    node["cells"] = ncell;
    node["cell_length"] = cell_length_;
    node["timestep"] = timestep_;
}

/**
 * initialize event list
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::init_event_list()
{
    step_ = 0;

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
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::compute_collision_event(const unsigned int n, cell_type const& cell)
{
    double dt = std::numeric_limits<double>::max();
    int n2 = -1;

    // iterate over particles in cell
    BOOST_FOREACH(unsigned int j, cell) {
	// skip same particle if in same cell
	if (j == n)
	    continue;

	// particle distance vector at time of first particle
	vector_type dr = part[j].r + part[j].v * (part[n].t - part[j].t) - part[n].r;
	// enforce periodic boundary conditions
	dr -= round(dr / box_) * box_;
	// velocity difference at given time
	vector_type dv = part[j].v - part[n].v;

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
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::compute_cell_event(const unsigned int n)
{
    vector_type dt3(std::numeric_limits<double>::max());
    double dt = std::numeric_limits<double>::max();
    cell_index cell2 = part[n].cell;

    for (unsigned int d = 0; d < dimension; ++d) {
	if (part[n].v[d] < 0.) {
	    dt3[d] = (part[n].cell[d] * cell_length_ - part[n].r[d]) / part[n].v[d];
	    cell2[d] = (cell2[d] + ncell - 1) % ncell;
	}
	else if (part[n].v[d] > 0.) {
	    dt3[d] = ((part[n].cell[d] + 1) * cell_length_ - part[n].r[d]) / part[n].v[d];
	    cell2[d] = (cell2[d] + 1) % ncell;
	}
	dt = std::min(dt, dt3[d]);
    }

    if (dt < event_list[n].t - part[n].t) {
	// generate cell boundary event
	event_list[n].t = part[n].t + dt;
	event_list[n].type = event::CELL;
	for (unsigned int d = 0; d < dimension; ++d) {
	    event_list[n].cell2[d] = (dt3[d] == dt) ? cell2[d] : part[n].cell[d];
	}
    }
}

/**
 * schedule next particle event starting at given time
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::schedule_event(const unsigned int n)
{
    // upper boundary for time of next particle event
    event_list[n].t = std::numeric_limits<double>::max();

    // compute next cell boundary event
    compute_cell_event(n);
    // compute next collision event with particles of this cell
    compute_collision_event(n, cell_(part[n].cell));

    // compute next collision event with particles of neighbour cells
    cell_index i;
    for (i[0] = -1; i[0] <= 1; ++i[0]) {
	for (i[1] = -1; i[1] <= 1; ++i[1]) {
	    if (dimension == 3) {
		for (i[2] = -1; i[2] <= 1; ++i[2]) {
		    cell_index j;
		    for (int d = 0; d < dimension; ++d) {
			j[d] = (part[n].cell[d] + ncell + i[d]) % ncell;
		    }
		    compute_collision_event(n, cell_(j));
		}
	    }
	    else {
		cell_index j;
		for (int d = 0; d < dimension; ++d) {
		    j[d] = (part[n].cell[d] + ncell + i[d]) % ncell;
		}
		compute_collision_event(n, cell_(j));
	    }
	}
    }

    // schedule particle event
    event_queue.push(event_queue_item(event_list[n].t, n));
}

/*
 * process particle collision event
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::process_collision_event(const unsigned int n1)
{
    const vector_type dr1 = part[n1].v * (event_list[n1].t - part[n1].t);
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

    const vector_type dr2 = part[n2].v * (event_list[n1].t - part[n2].t);
    // update periodically extended particle position
    part[n2].R += dr2;
    // update periodically reduced particle position to given time
    part[n2].r += dr2;
    // update particle time
    part[n2].t = event_list[n1].t;

    // particle distance vector
    vector_type dr = part[n2].r - part[n1].r;
    // enforce periodic boundary conditions
    dr -= round(dr / box_) * box_;
    // velocity difference before collision
    vector_type dv = part[n1].v - part[n2].v;
    // velocity difference after collision without dissipation
    dv = dr * (dr * dv) / (dr * dr);

    // update velocities to current simulation time
    part[n1].v -= dv;
    part[n2].v += dv;

    // add contribution to impulsive limit of the virial expression sum
    virial += dr * dv;

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
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::process_cell_event(const unsigned int n)
{
    const vector_type dr = part[n].v * (event_list[n].t - part[n].t);
    // update periodically extended particle position
    part[n].R += dr;
    // update periodically reduced particle position to given time
    part[n].r += dr;
    // enforce periodic boundary conditions
    for (unsigned int d = 0; d < dimension; ++d) {
	if (part[n].cell[d] == ncell - 1 && event_list[n].cell2[d] == 0)
	    part[n].r[d] -= box_;
	if (part[n].cell[d] == 0 && event_list[n].cell2[d] == ncell - 1)
	    part[n].r[d] += box_;
    }
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
template <int dimension>
typename hardsphere<hardsphere_impl<dimension> >::cell_index
hardsphere<hardsphere_impl<dimension> >::compute_cell(vector_type const& r)
{
    cell_index cell;
    for (int i = 0; i < dimension; ++i) {
	cell[i] = (int)(r[i] / cell_length_) % ncell;
    }
    return cell;
}

/**
 * advance phase space state to given sample time
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::mdstep()
{
    // nanosecond resolution process times
    boost::array<timespec, 2> t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[0]);

    // impulsive limit of the virial expression sum
    virial = 0.;

    // advance simulation time
    const double sample_time = ++step_ * timestep_;
    // process particle event queue till sample time
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
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[1]);

    virial /= npart;

    // CPU ticks for event queue processing
    m_times["event_queue"] += t[1] - t[0];
}

/**
 * advance phase space state to given sample time
 */
template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::sample(host_sample_type& sample) const
{
    typedef host_sample_type sample_type;
    typedef typename sample_type::position_sample_vector position_sample_vector;
    typedef typename sample_type::position_sample_ptr position_sample_ptr;
    typedef typename sample_type::velocity_sample_vector velocity_sample_vector;
    typedef typename sample_type::velocity_sample_ptr velocity_sample_ptr;

    // nanosecond resolution process times
    boost::array<timespec, 2> t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[0]);

    // allocate memory for phase space sample
    position_sample_ptr r(new position_sample_vector);
    velocity_sample_ptr v(new velocity_sample_vector);
    sample.r.push_back(r);
    sample.v.push_back(v);
    r->reserve(npart);
    v->reserve(npart);

    double const sample_time = step_ * static_cast<double>(timestep_);
    typename std::vector<particle>::const_iterator p;

    for (p = part.begin(); p != part.end(); ++p) {
	// periodically extended particle position
	r->push_back(p->R + p->v * (sample_time - p->t));
	// particle velocity
	v->push_back(p->v);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[1]);

    // CPU ticks for phase space sampling
    m_times["sample"] += t[1] - t[0];
}

template <int dimension>
void hardsphere<hardsphere_impl<dimension> >::sample(energy_sample_type& sample) const
{
    typename std::vector<particle>::const_iterator p;

    // mean potential energy per particle
    sample.en_pot = 0;
    // mean virial equation sum per particle
    sample.virial = virial;

    for (p = part.begin(), sample.vv = 0, sample.v_cm = 0; p != part.end(); ++p) {
	sample.vv += p->v * p->v;
	sample.v_cm += p->v;
    }
    // mean squared velocity per particle
    sample.vv /= npart;
    // mean velocity per particle
    sample.v_cm /= npart;
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_HARDSPHERE_HPP */
