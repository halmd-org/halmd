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
    typedef typename sample_type::position_vector position_vector;
    typedef typename sample_type::velocity_vector velocity_vector;

    /**
     * MD simulation particle
     */
    struct particle
    {
	typedef boost::reference_wrapper<particle> ref;

	/** particle position */
	vector_type r;
	/** particle velocity */
	vector_type v;
	/** particle force */
	vector_type f;
	/** particle number */
	unsigned int tag;
	/** particle neighbours list */
	std::vector<ref> neighbour;
    };

    typedef typename std::vector<typename particle::ref> cell_list;
    typedef typename cell_list::iterator cell_list_iterator;
    typedef typename cell_list::const_iterator cell_list_const_iterator;
    typedef boost::array<int, dimension> cell_index;
    typedef boost::multi_array<cell_list, dimension> cell_lists;

public:
    /** set number of particles */
    void particles(unsigned int value);
    /** set neighbour list skin */
    void nbl_skin(float value);

    /** set system state from phase space sample */
    void sample(sample_visitor visitor);
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
    /** returns trajectory sample */
    sample_type const& sample() const { return m_sample; }
    /** returns number of cells per dimension */
    int cells() const { return ncell; }
    /** returns cell length */
    double cell_length() const { return cell_length_; }

    /** MD simulation step */
    void mdstep();

    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

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
    /** compute C²-smooth potential */
    void compute_smooth_potential(double r, double fval, double pot);
    /** first leapfrog step of integration of equations of motion */
    void leapfrog_half();
    /** second leapfrog step of integration of equations of motion */
    void leapfrog_full();
    /** random collision with heat bath */
    void anderson_thermostat();

private:
    using _Base::npart;
    using _Base::density_;
    using _Base::box_;
    using _Base::timestep_;
    using _Base::r_cut;
    using _Base::rr_cut;
    using _Base::en_cut;
    using _Base::r_smooth;
    using _Base::rri_smooth;
    using _Base::thermostat_nu;
    using _Base::thermostat_temp;

    using _Base::m_sample;
    using _Base::m_times;

    /** particles */
    std::vector<particle> part;
    /** cell lists */
    cell_lists cell;
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

/**
 * set number of particles in system
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::particles(unsigned int value)
{
    _Base::particles(value);

    try {
	m_sample.r.resize(npart);
	m_sample.v.resize(npart);
	part.resize(npart);
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate phase space state");
    }
}

/**
 * set system state from phase space sample
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::sample(sample_visitor visitor)
{
    _Base::sample(visitor);

    for (unsigned int i = 0; i < npart; ++i) {
	part[i].r = m_sample.r[i];
	part[i].v = m_sample.v[i];
	part[i].tag = i;
    }
    // update cell lists
    update_cells();
    // update Verlet neighbour lists
    update_neighbours();
    // reset sum over maximum velocity magnitudes to zero
    v_max_sum = 0.;
    // calculate forces, potential energy and virial equation sum
    compute_forces();
}

template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::nbl_skin(float value)
{
    r_skin = value;
    LOG("neighbour list skin: " << r_skin);

    // cutoff radius with neighbour list skin
    r_cut_skin = r_skin + r_cut;
    // squared cutoff radius with neighbour list skin
    rr_cut_skin = r_cut_skin * r_cut_skin;

    // number of cells per dimension
    ncell = std::floor(box_ / r_cut_skin);
    LOG("number of cells per dimension: " << ncell);

    if (ncell < 3) {
	throw exception("requires at least 3 cells per dimension");
    }

    // create empty cell lists
    cell_index size;
    std::fill(size.begin(), size.end(), ncell);
    cell.resize(size);

    // derive cell length from integer number of cells per dimension
    cell_length_ = box_ / ncell;
    LOG("cell length: " << cell_length_);
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
	part[i].tag = i;
	vector_type& r = part[i].r;
	// compose primitive vectors from 1-dimensional index
	if (dimension == 3) {
	    r[0] = ((i >> 2) % n) + ((i ^ (i >> 1)) & 1) / 2.;
	    r[1] = ((i >> 2) / n % n) + (i & 1) / 2.;
	    r[2] = ((i >> 2) / n / n) + (i & 2) / 4.;
	}
	else {
	    r[0] = ((i >> 1) % n) + (i & 1) / 2.;
	    r[1] = ((i >> 1) / n) + (i & 1) / 2.;
	}
	r *= a;
    }

    // update cell lists
    update_cells();
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
    vector_type v_cm = 0;
    // maximum squared velocity
    double vv_max = 0;

    BOOST_FOREACH(particle& p, part) {
	// generate random Maxwell-Boltzmann distributed velocity
	rng_.gaussian(p.v, value);
	v_cm += p.v;
	// initialize force to zero for first leapfrog half step
	p.f = 0;
    }

    v_cm /= npart;

    BOOST_FOREACH(particle& p, part) {
	// set center of mass velocity to zero
	p.v -= v_cm;

	vv_max = std::max(vv_max, p.v * p.v);
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
    // create empty cell lists
    cell_index size;
    std::fill(size.begin(), size.end(), ncell);
    cell = cell_lists(size);
    // add particles to cells
    BOOST_FOREACH(particle& p, part) {
	compute_cell(p.r).push_back(boost::ref(p));
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
	if (same_cell && p2.tag <= p1.tag)
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
    BOOST_FOREACH(particle& p, part) {
	p.f = 0;
    }

    // potential energy
    m_sample.en_pot = 0;
    // virial equation sum
    m_sample.virial = 0;

    BOOST_FOREACH(particle& p1, part) {
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
	    double rri = 1 / rr;
	    double r6i = rri * rri * rri;
	    double fval = 48 * rri * r6i * (r6i - 0.5);
	    double pot = 4 * r6i * (r6i - 1) - en_cut;

	    if (r_smooth > 0) {
		compute_smooth_potential(std::sqrt(rr), fval, pot);
	    }

	    // add force contribution to both particles
	    p1.f += r * fval;
	    p2.f -= r * fval;

	    // add contribution to potential energy
	    m_sample.en_pot += pot;
	    // add contribution to virial equation sum
	    m_sample.virial += rr * fval;
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
 * compute C²-smooth potential
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::compute_smooth_potential(double r, double fval, double pot)
{
    double y = r - r_cut;
    double x2 = y * y * rri_smooth;
    double x4 = x2 * x2;
    double x4i = 1 / (1 + x4);
    // smoothing function
    double h0_r = x4 * x4i;
    // first derivative of smoothing function
    double h1_r = 4 * y * x2 * x4i * x4i;
    // apply smoothing function to obtain C¹ force function
    fval = h0_r * fval - h1_r * pot / r;
    // apply smoothing function to obtain C² potential function
    pot = h0_r * pot;
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::leapfrog_half()
{
    BOOST_FOREACH(particle& p, part) {
	// half step velocity
	p.v += p.f * (timestep_ / 2.);
	// full step position
	p.r += p.v * timestep_;
	// copy to phase space sample
	m_sample.r[p.tag] = p.r;
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

    BOOST_FOREACH(particle& p, part) {
	// full step velocity
	p.v += p.f * (timestep_ / 2.);
	// copy to phase space sample
	m_sample.v[p.tag] = p.v;

	vv_max = std::max(vv_max, p.v * p.v);
    }

    // add to sum over maximum velocity magnitudes since last neighbour lists update
    v_max_sum += std::sqrt(vv_max);
}

/**
 * random collisions with heat bath
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::anderson_thermostat()
{
    BOOST_FOREACH(particle& p, part) {
	if (rng_.uniform() < (thermostat_nu * timestep_)) {
	    rng_.gaussian(p.v, thermostat_temp);
	}
    }
}

/**
 * MD simulation step
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::mdstep()
{
    // nanosecond resolution process times
    boost::array<timespec, 7> t;

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
    if (thermostat_nu > 0) {
	anderson_thermostat();
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[5]);

	m_times["anderson_thermostat"] += t[5] - t[4];
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[6]);

    m_times["update_forces"] += t[3] - t[2];
    m_times["velocity_verlet"] += (t[1] - t[0]) + (t[4] - t[3]);
    m_times["mdstep"] += t[6] - t[0];
}

/**
 * write parameters to HDF5 parameter group
 */
template <int dimension>
void ljfluid<ljfluid_impl_host<dimension> >::param(H5param& param) const
{
    _Base::param(param);

    H5xx::group node(param["mdsim"]);
    node["cells"] = ncell;
    node["cell_length"] = cell_length_;
    node["neighbour_skin"] = r_skin;
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_HOST_HPP */
