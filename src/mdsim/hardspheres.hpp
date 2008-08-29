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
#include <boost/function.hpp>
#include <boost/multi_array.hpp>
#include <list>
#include <queue>
#include <vector>
#include "config.hpp"
#include "gsl_rng.hpp"
#include "perf.hpp"

namespace mdsim
{

/**
 * Hard Spheres simulation
 */
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
	hvector r;
	/** periodically extended particle position */
	hvector R;
	/** particle velocity */
	hvector v;
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

    /** trajectory sample visitor type */
    typedef boost::function<void (std::vector<hvector>&, std::vector<hvector>&)> trajectory_sample_visitor;
    /** MD simulation sample visitor type */
    typedef boost::function<void (std::vector<hvector> const&, std::vector<hvector> const&, std::vector<hvector> const&, double)> mdsim_sample_visitor;

public:
    hardspheres() : step_(0) {}
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
    /** set simulation timestep */
    void timestep(double value);

    /** set system state from phase space sample */
    void restore(trajectory_sample_visitor visitor);
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
    unsigned int const& particles() const { return npart; }
    /** returns pair separation at which particle collision occurs */
    double const& pair_separation() const { return pair_sep_; }
    /** returns particle density */
    double const& density() const { return density_; }
    /** returns periodic box length */
    double const& box() const { return box_; }
    /** returns number of cells per dimension */
    unsigned int const& cells() const { return ncell; }
    /** returns cell length */
    double const& cell_length() { return cell_length_; }
    /** returns simulation timestep */
    double const& timestep() { return timestep_; }
    /** returns and resets CPU tick statistics */
    perf_counters times();

    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;

    /** advance phase space state to given sample time */
    void mdstep();
    /** sample trajectory */
    void sample(mdsim_sample_visitor visitor) const;

private:
    /** schedule next particle event starting at given time */
    void schedule_event(unsigned int n);
    /** process particle collision event */
    void process_collision_event(unsigned int n);
    /** process cell boundary event */
    void process_cell_event(unsigned int n);
    /** returns cell which a particle belongs to */
    cell_index compute_cell(hvector const& r);
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

    /** periodically reduced particle positions at sample time */
    std::vector<hvector> r_;
    /** periodically extended particle positions at sample time */
    std::vector<hvector> R_;
    /** particle velocities at sample time */
    std::vector<hvector> v_;
    /** impulsive limit of the virial expression sum */
    double virial_;

    /** random number generator */
    rng::gsl::gfsr4 rng_;
    /** squared pair separation */
    double pair_sep_sq;

    /** CPU tick statistics */
    perf_counters m_times;
};

} // namespace mdsim

#endif /* ! MDSIM_HARDSPHERES_HPP */
