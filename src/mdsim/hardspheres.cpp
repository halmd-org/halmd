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

#include <algorithm>
#include <boost/foreach.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <sys/times.h>
#include "H5xx.hpp"
#include "exception.hpp"
#include "hardspheres.hpp"
#include "log.hpp"
#include "timer.hpp"

#define foreach BOOST_FOREACH

namespace mdsim
{

/**
 * set number of particles in system
 */
void hardspheres::particles(unsigned int value)
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
void hardspheres::pair_separation(double value)
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
void hardspheres::density(double value)
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
void hardspheres::box(double value)
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
void hardspheres::init_cell()
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
 * set simulation timestep
 */
void hardspheres::timestep(double value)
{
    timestep_ = value;
    LOG("simulation timestep: " << timestep_);
}

/**
 * set system state from phase space sample
 */
void hardspheres::restore(trajectory_sample_visitor visitor)
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
void hardspheres::rng(unsigned int seed)
{
    rng_.set(seed);
    LOG("initializing random number generator with seed: " << seed);
}

/**
 * initialize random number generator from state
 */
void hardspheres::rng(rng::gsl::gfsr4::state_type const& state)
{
    rng_.restore(state);
    LOG("restoring random number generator from state");
}

/**
 * place particles on a face-centered cubic (fcc) lattice
 */
void hardspheres::lattice()
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
#ifdef DIM_3D
	// compose primitive vectors from 1-dimensional index
	part[i].r[0] = ((i >> 2) % n) + ((i ^ (i >> 1)) & 1) / 2.;
	part[i].r[1] = ((i >> 2) / n % n) + (i & 1) / 2.;
	part[i].r[2] = ((i >> 2) / n / n) + (i & 2) / 4.;
#else
	// compart[i]ose primitive vectors from 1-dimensional index
	part[i].r[0] = ((i >> 1) % n) + (i & 1) / 2.;
	part[i].r[1] = ((i >> 1) / n) + (i & 1) / 2.;
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
void hardspheres::temperature(double value)
{
    LOG("initializing velocities from Maxwell-Boltzmann distribution at temperature: " << value);

    // center of mass velocity
    hvector v_cm = 0.;

    foreach (particle& p, part) {
	// generate random Maxwell-Boltzmann distributed velocity
	rng_.gaussian(p.v[0], p.v[1], value);
#ifdef DIM_3D
	// Box-Muller transformation strictly generates 2 variates at once
	rng_.gaussian(p.v[1], p.v[2], value);
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
 * write parameters to HDF5 parameter group
 */
void hardspheres::attrs(H5::Group const& param) const
{
    H5xx::group node(param.createGroup("mdsim"));
    node["dimension"] = dimension;
    node["particles"] = npart;
    node["pair_separation"] = pair_sep_;
    node["cells"] = ncell;
    node["cell_length"] = cell_length_;
    node["density"] = density_;
    node["box_length"] = box_;
    node["timestep"] = timestep_;
}

/**
 * initialize event list
 */
void hardspheres::init_event_list()
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
void hardspheres::compute_collision_event(const unsigned int n, cell_type const& cell)
{
    double dt = std::numeric_limits<double>::max();
    int n2 = -1;

    // iterate over particles in cell
    foreach (unsigned int j, cell) {
	// skip same particle if in same cell
	if (j == n)
	    continue;

	// particle distance vector at time of first particle
	hvector dr = part[j].r + part[j].v * (part[n].t - part[j].t) - part[n].r;
	// enforce periodic boundary conditions
	dr -= round(dr / box_) * box_;
	// velocity difference at given time
	hvector dv = part[j].v - part[n].v;

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
void hardspheres::compute_cell_event(const unsigned int n)
{
    hvector dt3(std::numeric_limits<double>::max());
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
void hardspheres::schedule_event(const unsigned int n)
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
void hardspheres::process_collision_event(const unsigned int n1)
{
    const hvector dr1 = part[n1].v * (event_list[n1].t - part[n1].t);
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

    const hvector dr2 = part[n2].v * (event_list[n1].t - part[n2].t);
    // update periodically extended particle position
    part[n2].R += dr2;
    // update periodically reduced particle position to given time
    part[n2].r += dr2;
    // update particle time
    part[n2].t = event_list[n1].t;

    // particle distance vector
    hvector dr = part[n2].r - part[n1].r;
    // enforce periodic boundary conditions
    dr -= round(dr / box_) * box_;
    // velocity difference before collision
    hvector dv = part[n1].v - part[n2].v;
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
void hardspheres::process_cell_event(const unsigned int n)
{
    const hvector dr = part[n].v * (event_list[n].t - part[n].t);
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
hardspheres::cell_index hardspheres::compute_cell(hvector const& r)
{
    hvector cellf = r / cell_length_;
#ifdef DIM_3D
    cell_index cell = {{ int(cellf[0]), int(cellf[1]), int(cellf[2]) }};
#else
    cell_index cell = {{ int(cellf[0]), int(cellf[1]) }};
#endif
    return cell;
}

/**
 * advance phase space state to given sample time
 */
void hardspheres::mdstep()
{
    // nanosecond resolution process times
    boost::array<timespec, 5> t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[0]);

    // impulsive limit of the virial expression sum
    virial_ = 0.;

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

    virial_ /= npart;

    // sample phase space at given time
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[1]);
    for (unsigned int i = 0; i < npart; ++i) {
	const hvector dr = part[i].v * (sample_time - part[i].t);
	// periodically extended particle position
	R_[i] = part[i].R + dr;
	// periodically reduced particle position
	r_[i] = part[i].r + dr;
	// enforce periodic boundary conditions
	r_[i] -= floor(r_[i] / box_) * box_;
	// particle velocity
	v_[i] = part[i].v;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[2]);

    // CPU ticks for MD simulation step
    m_times[0] += t[2] - t[0];
    // CPU ticks for event queue processing
    m_times[1] += t[1] - t[0];
    // CPU ticks for phase space sampling
    m_times[2] += t[2] - t[1];
}

/**
 * sample trajectory
 */
void hardspheres::sample(mdsim_sample_visitor visitor) const
{
    visitor(r_, R_, v_, virial_);
}
 
/*
 * returns and resets CPU tick statistics
 */
perf_counters hardspheres::times()
{
    perf_counters times(m_times);
    // reset performance counters
    for (unsigned int i = 0; i < m_times.size(); ++i) {
	m_times[i].clear();
    }
    return times;
}

} // namespace mdsim
