/* Lennard-Jones fluid simulation
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

#ifndef LJGPU_MDSIM_LJFLUID_HOST_HPP
#define LJGPU_MDSIM_LJFLUID_HOST_HPP

#include <algorithm>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/ref.hpp>
#include <cmath>
#include <iostream>
#include <list>
#include <ljgpu/mdsim/ljfluid_base.hpp>
#include <ljgpu/rng/gsl_rng.hpp>
#include <ljgpu/util/timer.hpp>
#include <sys/times.h>
#include <vector>

namespace ljgpu
{

template <typename ljfluid_impl>
class ljfluid;

template<int dimension>
class ljfluid<ljfluid_impl_host<dimension> >
    : public ljfluid_base<ljfluid_impl_host<dimension> >
{
public:
    typedef ljfluid_base<ljfluid_impl_host<dimension> > _Base;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::sample_type sample_type;
    typedef typename sample_type::sample_visitor sample_visitor;

    /**
     * MD simulation particle
     */
    struct particle
    {
	/** particle position */
	vector_type r;
	/** particle velocity */
	vector_type v;
	/** particle number tag */
	int n;

	/** particle force */
	vector_type f;
	/** particle neighbours list */
	std::vector<boost::reference_wrapper<particle> > neighbour;

	particle(vector_type const& r, int n) : r(r), n(n) {}
	particle(vector_type const& r, vector_type const& v, int n) : r(r), v(v), n(n) {}
	particle() {}
    };

    typedef typename std::list<particle> cell_list;
    typedef typename cell_list::iterator cell_list_iterator;
    typedef typename cell_list::const_iterator cell_list_const_iterator;
    typedef boost::array<int, dimension> cell_index;

public:
    /** initialise fluid from program options */
    ljfluid(options const& opt);

    using _Base::density;
    using _Base::box;
    using _Base::timestep;
#ifdef USE_POTENTIAL_SMOOTHING
    using _Base::potential_smoothing;
#endif

    /** set number of particles */
    void particles(unsigned int value);
    /** set potential cutoff radius */
    void cutoff_radius(float_type value);
    /** initialize cell lists */
    void init_cell();

    /** set system state from phase space sample */
    void restore(sample_visitor visitor);
    /** initialize random number generator with seed */
    void rng(unsigned int seed);
    /** initialize random number generator from state */
    void rng(gsl::gfsr4::state_type const& state);
    /** place particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(double value);

    /** returns number of particles */
    unsigned int particles() const { return npart; }
    /** returns potential cutoff radius */
    float_type cutoff_radius() const { return r_cut; }
    /** returns number of cells per dimension */
    int cells() const { return ncell; }
    /** returns cell length */
    double cell_length() const { return cell_length_; }

    /** MD simulation step */
    void mdstep();
    /** ljfluid GPU compat */
    void stream() {}
    /** ljfluid GPU compat */
    void copy() {}

    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;

private:
    /** update cell lists */
    void update_cells();
    /** returns cell list which a particle belongs to */
    cell_list& compute_cell(vector_type r);
    /** update neighbour lists */
    void update_neighbours();
    /** update neighbour lists for a single cell */
    void update_cell_neighbours(cell_index const& i);
    /** update neighbour list of particle */
    template <bool same_cell> void compute_cell_neighbours(particle& p, cell_list& c);
    /** compute Lennard-Jones forces */
    void compute_forces();
    /** first leapfrog step of integration of equations of motion */
    void leapfrog_half();
    /** second leapfrog step of integration of equations of motion */
    void leapfrog_full();

private:
    using _Base::npart;
    using _Base::density_;
    using _Base::box_;
    using _Base::timestep_;
    using _Base::r_cut;
    using _Base::rr_cut;
    using _Base::en_cut;
#ifdef USE_POTENTIAL_SMOOTHING
    using _Base::r_smooth;
    using _Base::rri_smooth;
#endif
    using _Base::m_sample;
    using _Base::m_times;

    /** cell lists */
    boost::multi_array<cell_list, dimension> cell;
    /** random number generator */
    gsl::gfsr4 rng_;

    /** number of cells per dimension */
    int ncell;
    /** cell length */
    float_type cell_length_;
    /** neighbour list skin */
    float_type r_skin;
    /** cutoff radius with neighbour list skin */
    float_type r_cut_skin;
    /** squared cutoff radius with neighbour list skin */
    float_type rr_cut_skin;
    /** sum over maximum velocity magnitudes since last neighbour lists update */
    float_type v_max_sum;
};

template <int dimension>
ljfluid<ljfluid_impl_host<dimension> >::ljfluid(options const& opt)
{
    LOG("positional coordinates dimension: " << dimension);

    particles(opt["particles"].as<unsigned int>());
    if (opt["density"].defaulted() && !opt["box-length"].empty()) {
	box(opt["box-length"].as<float>());
    }
    else {
	density(opt["density"].as<float>());
    }
    cutoff_radius(opt["cutoff"].as<float>());
#ifdef USE_POTENTIAL_SMOOTHING
    potential_smoothing(opt["smoothing"].as<float>());
#endif
    timestep(opt["timestep"].as<float>());
    init_cell();
}

/**
 * set potential cutoff radius
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::cutoff_radius(float_type value)
{
    _Base::cutoff_radius(value);

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
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::particles(unsigned int value)
{
    _Base::particles(value);

    try {
	m_sample.R.resize(npart);
	m_sample.v.resize(npart);
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate phase space state");
    }
}

/**
 * set system state from phase space sample
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::restore(sample_visitor visitor)
{
    // set system state from phase space sample
    visitor(m_sample.R, m_sample.v);

    for (unsigned int i = 0; i < npart; ++i) {
	// add particle to appropriate cell list
	compute_cell(m_sample.R[i]).push_back(particle(m_sample.R[i], m_sample.v[i], i));
    }

    // update Verlet neighbour lists
    update_neighbours();
    // reset sum over maximum velocity magnitudes to zero
    v_max_sum = 0.;
    // calculate forces, potential energy and virial equation sum
    compute_forces();
}

/**
 * initialize cell lists
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::init_cell()
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
	cell_index size;
	std::fill(size.begin(), size.end(), ncell);
	cell.resize(size);
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate initial cell lists");
    }
}

/**
 * initialize random number generator with seed
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::rng(unsigned int seed)
{
    rng_.set(seed);
    LOG("initializing random number generator with seed: " << seed);
}

/**
 * initialize random number generator from state
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::rng(gsl::gfsr4::state_type const& state)
{
    rng_.restore(state);
    LOG("restoring random number generator from state");
}

/**
 * place particles on a face-centered cubic (fcc) lattice
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::lattice()
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
	vector_type r(a);
	// compose primitive vectors from 1-dimensional index
	if (dimension == 3) {
	    r[0] *= ((i >> 2) % n) + ((i ^ (i >> 1)) & 1) / 2.;
	    r[1] *= ((i >> 2) / n % n) + (i & 1) / 2.;
	    r[2] *= ((i >> 2) / n / n) + (i & 2) / 4.;
	}
	else {
	    r[0] *= ((i >> 1) % n) + (i & 1) / 2.;
	    r[1] *= ((i >> 1) / n) + (i & 1) / 2.;
	}
	// add particle to appropriate cell list
	compute_cell(r).push_back(particle(r, i));
	// copy position to sorted particle list
	m_sample.R[i] = r;
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
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::temperature(double value)
{
    LOG("initializing velocities from Maxwell-Boltzmann distribution at temperature: " << value);

    // center of mass velocity
    vector_type v_cm = 0.;
    // maximum squared velocity
    double vv_max = 0.;

    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	BOOST_FOREACH(particle& p, *it) {
	    // generate random Maxwell-Boltzmann distributed velocity
	    rng_.gaussian(p.v[0], p.v[1], value);
	    if (dimension == 3) {
		// Box-Muller transformation strictly generates 2 variates at once
		rng_.gaussian(p.v[1], p.v[2], value);
	    }
	    v_cm += p.v;

	    // initialize force to zero for first leapfrog half step
	    p.f = 0.;
	}
    }

    v_cm /= npart;

    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	BOOST_FOREACH(particle& p, *it) {
	    // set center of mass velocity to zero
	    p.v -= v_cm;
	    // copy velocity to sorted particle list
	    m_sample.v[p.n] = p.v;

	    vv_max = std::max(vv_max, p.v * p.v);
	}
    }

    // initialize sum over maximum velocity magnitudes since last neighbour lists update
    v_max_sum = std::sqrt(vv_max);
}

/**
 * update cell lists
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::update_cells()
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
template <int dimension>
typename ljfluid<ljfluid_impl_host<dimension> >::cell_list&
ljfluid<ljfluid_impl_host<dimension> >::compute_cell(vector_type r)
{
    r = (r - floor(r / box_) * box_) / cell_length_;
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
    cell_index idx;
    for (int i = 0; i < dimension; ++i) {
	idx[i] = (unsigned int)(r[i]) % ncell;
    }
    return cell(idx);
}

/**
 * update neighbour lists
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::update_neighbours()
{
    cell_index i;
    for (i[0] = 0; i[0] < ncell; ++i[0]) {
	for (i[1] = 0; i[1] < ncell; ++i[1]) {
	    if (dimension == 3) {
		for (i[2] = 0; i[2] < ncell; ++i[2]) {
		    update_cell_neighbours(i);
		}
	    }
	    else {
		update_cell_neighbours(i);
	    }
	}
    }
}

/**
 * update neighbour lists for a single cell
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::update_cell_neighbours(cell_index const& i)
{
    BOOST_FOREACH(particle& p, cell(i)) {
	// empty neighbour list of particle
	p.neighbour.clear();

	cell_index j;
	for (j[0] = -1; j[0] <= 1; ++j[0]) {
	    for (j[1] = -1; j[1] <= 1; ++j[1]) {
		if (dimension == 3) {
		    for (j[2] = -1; j[2] <= 1; ++j[2]) {
			// visit half of 26 neighbour cells due to pair potential
			if (j[0] == 0 && j[1] == 0 && j[2] == 0) {
			    goto out;
			}
			// update neighbour list of particle
			cell_index k;
			for (int n = 0; n < dimension; ++n) {
			    k[n] = (i[n] + ncell + j[n]) % ncell;
			}
			compute_cell_neighbours<false>(p, cell(k));
		    }
		}
		else {
		    // visit half of 8 neighbour cells due to pair potential
		    if (j[0] == 0 && j[1] == 0) {
			goto out;
		    }
		    // update neighbour list of particle
		    boost::array<unsigned int, dimension> k;
		    for (int n = 0; n < dimension; ++n) {
			k[n] = (i[n] + ncell + j[n]) % ncell;
		    }
		    compute_cell_neighbours<false>(p, cell(k));
		}
	    }
	}
out:
	// visit this cell
	compute_cell_neighbours<true>(p, cell(i));
    }
}

/**
 * update neighbour list of particle
 */
template <int dimension>
template <bool same_cell>
void ljfluid<ljfluid_impl_host<dimension> >::compute_cell_neighbours(particle& p1, cell_list& c)
{
    BOOST_FOREACH(particle& p2, c) {
	// skip identical particle and particle pair permutations if same cell
	if (same_cell && p2.n <= p1.n)
	    continue;

	// particle distance vector
	vector_type r = p1.r - p2.r;
	// enforce periodic boundary conditions
	r -= round(r / box_) * box_;
	// squared particle distance
	double rr = r * r;

	// enforce cutoff radius with neighbour list skin
	if (rr >= rr_cut_skin)
	    continue;

	// add particle to neighbour list
	p1.neighbour.push_back(boost::ref(p2));
    }
}

/**
 * compute Lennard-Jones forces
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::compute_forces()
{
    // initialize particle forces to zero
    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	BOOST_FOREACH(particle& p, *it) {
	    p.f = 0.;
	}
    }

    // potential energy
    m_sample.en_pot = 0.;
    // virial equation sum
    m_sample.virial = 0.;

    // iterate over all particles
    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	BOOST_FOREACH(particle& p1, *it) {
	    // calculate pairwise Lennard-Jones force with neighbour particles
	    BOOST_FOREACH(particle& p2, p1.neighbour) {
		// particle distance vector
		vector_type r = p1.r - p2.r;
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
		p1.f += r * fval;
		p2.f -= r * fval;

		// add contribution to potential energy
		m_sample.en_pot += 4. * r6i * (r6i - 1.) - en_cut;
		// add contribution to virial equation sum
		m_sample.virial += rr * fval;
	    }
	}
    }

    m_sample.en_pot /= npart;
    m_sample.virial /= npart;

    // ensure that system is still in valid state
    if (std::isinf(m_sample.en_pot)) {
	throw exception("potential energy diverged due to excessive timestep or density");
    }
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::leapfrog_half()
{
    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	BOOST_FOREACH(particle& p, *it) {
	    // half step velocity
	    p.v += p.f * (timestep_ / 2.);
	    // full step position
	    p.r += p.v * timestep_;
	    // copy position to sorted particle list
	    m_sample.R[p.n] = p.r;
	}
    }
}

/**
 * second leapfrog step of integration of equations of motion
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::leapfrog_full()
{
    // maximum squared velocity
    double vv_max = 0.;

    for (cell_list* it = cell.data(); it != cell.data() + cell.num_elements(); ++it) {
	BOOST_FOREACH(particle& p, *it) {
	    // full step velocity
	    p.v += p.f * (timestep_ / 2.);
	    // copy velocity to sorted particle list
	    m_sample.v[p.n] = p.v;

	    vv_max = std::max(vv_max, p.v * p.v);
	}
    }

    // add to sum over maximum velocity magnitudes since last neighbour lists update
    v_max_sum += std::sqrt(vv_max);
}

/**
 * MD simulation step
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::mdstep()
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

	m_times["update_cells"] += t[2] - t[1];
	m_times["update_neighbours"] += t[3] - t[2];
    }

    // calculate forces, potential energy and virial equation sum
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[2]);
    compute_forces();
    // calculate velocities
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[3]);
    leapfrog_full();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[4]);

    m_times["update_forces"] += t[3] - t[2];
    m_times["velocity_verlet"] += (t[1] - t[0]) + (t[4] - t[3]);
    m_times["mdstep"] += t[4] - t[0];
}

/**
 * write parameters to HDF5 parameter group
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::attrs(H5::Group const& param) const
{
    _Base::attrs(param);

    H5xx::group node(param.openGroup("mdsim"));
    node["cells"] = cells();
    node["cell_length"] = cell_length();
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_HOST_HPP */
