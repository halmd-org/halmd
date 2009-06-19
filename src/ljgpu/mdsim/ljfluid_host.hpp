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
#include <boost/bind.hpp>
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

#define foreach BOOST_FOREACH
#define range boost::make_iterator_range

namespace ljgpu
{

template <typename ljfluid_impl, int dimension>
class ljfluid;

template <int dimension>
class ljfluid<ljfluid_impl_host, dimension>
    : public ljfluid_base<ljfluid_impl_host, dimension>
{
public:
    typedef ljfluid_base<ljfluid_impl_host, dimension> _Base;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::host_sample_type host_sample_type;
    typedef typename _Base::energy_sample_type energy_sample_type;
    typedef typename _Base::virial_tensor virial_tensor;

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
	/** particle type */
	enum types { A = 0, B = 1 } type;
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
    template <typename T>
    void particles(T const& value);
    /** set neighbour list skin */
    void nbl_skin(float value);

    /** set system state from phase space sample */
    void state(host_sample_type& sample, float_type box);
    /** rescale particle velocities */
    void rescale_velocities(double coeff);
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
    /** returns number of cells per dimension */
    int cells() const { return ncell; }
    /** returns cell length */
    double cell_length() const { return cell_length_; }

    /** MD simulation step */
    void mdstep();
    /** sample phase space on host */
    void sample(host_sample_type& sample) const;
    /** sample thermodynamic equilibrium properties */
    void sample(energy_sample_type& sample) const;

    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

private:
    /** initialise velocities from Maxwell-Boltzmann distribution */
    void boltzmann(double temp);
    /** randomly assign particles types in a binary mixture */
    void random_binary_types();
    /** update cell lists */
    void update_cells();
    /** returns cell list which a particle belongs to */
    cell_list& compute_cell(vector_type r);
    /** update neighbour lists */
    template <bool binary>
    void update_neighbours();
    /** update neighbour lists for a single cell */
    template <bool binary>
    void update_cell_neighbours(cell_index const& i);
    /** update neighbour list of particle */
    template <bool same_cell, bool binary>
    void compute_cell_neighbours(particle& p, cell_list& c);
    /** compute Lennard-Jones forces */
    template <bool binary>
    void compute_forces();
    /** compute C²-smooth potential */
    template <bool binary>
    void compute_smooth_potential(float_type r, float_type& fval, float_type& pot, unsigned int type);
    /** first leapfrog step of integration of equations of motion */
    void leapfrog_half();
    /** second leapfrog step of integration of equations of motion */
    void leapfrog_full();

private:
    using _Base::npart;
    using _Base::mpart;
    using _Base::density_;
    using _Base::box_;
    using _Base::timestep_;
    using _Base::r_cut;
    using _Base::rr_cut;
    using _Base::en_cut;
    using _Base::sigma_;
    using _Base::sigma2_;
    using _Base::epsilon_;
    using _Base::r_smooth;
    using _Base::rri_smooth;
    using _Base::thermostat_steps;
    using _Base::thermostat_count;
    using _Base::thermostat_temp;

    using _Base::m_times;

    using _Base::mixture_;
    using _Base::potential_;

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
    /** cutoff radii with neighbour list skin */
    boost::array<float_type, 3> r_cut_skin;
    /** squared cutoff radii with neighbour list skin */
    boost::array<float_type, 3> rr_cut_skin;

    /** potential energy per particle */
    double en_pot;
    /** virial equation sum per particle */
    std::vector<virial_tensor> virial;
    /** sum over maximum velocity magnitudes since last neighbour lists update */
    float_type v_max_sum;
};

/**
 * set number of particles in system
 */
template <int dimension>
template <typename T>
void ljfluid<ljfluid_impl_host, dimension>::particles(T const& value)
{
    _Base::particles(value);

    try {
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
void ljfluid<ljfluid_impl_host, dimension>::state(host_sample_type& sample, float_type box)
{
    typedef typename host_sample_type::value_type sample_type;

    _Base::state(sample, box);

    using namespace boost::assign;
    boost::array<typename particle::types, 2> const types = list_of(particle::A)(particle::B);
    typename sample_type::position_sample_vector::const_iterator r;
    typename sample_type::velocity_sample_vector::const_iterator v;

    for (size_t i = 0, n = 0; n < npart; ++i) {
	for (r = sample[i].r->begin(), v = sample[i].v->begin(); r != sample[i].r->end(); ++r, ++v, ++n) {
	    part[n].type = types[i];
	    part[n].tag = n;
	    part[n].r = *r;
	    part[n].v = *v;
	}
    }

    // update cell lists
    update_cells();

    if (mixture_ == BINARY) {
	// update Verlet neighbour lists
	update_neighbours<true>();
	// calculate forces, potential energy and virial equation sum
	compute_forces<true>();
    }
    else {
	update_neighbours<false>();
	compute_forces<false>();
    }

    // reset sum over maximum velocity magnitudes to zero
    v_max_sum = 0;
}

template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::nbl_skin(float value)
{
    r_skin = value;
    LOG("neighbour list skin: " << r_skin);

    for (size_t i = 0; i < sigma_.size(); ++i) {
	r_cut_skin[i] = r_cut[i] + r_skin;
	rr_cut_skin[i] = std::pow(r_cut_skin[i], 2);
    }

    // number of cells per dimension
    ncell = static_cast<int>(box_ / *std::max_element(r_cut_skin.begin(), r_cut_skin.end()));
    LOG("number of cells per dimension: " << ncell);

    if (ncell < 3) {
	throw exception("less than least 3 cells per dimension");
    }

    // create empty cell lists
    cell_index size;
    std::fill(size.begin(), size.end(), ncell);
    cell.resize(size);

    // derive cell length from integer number of cells per dimension
    cell_length_ = box_ / ncell;
    LOG("cell length: " << cell_length_);
}

template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::rescale_velocities(double coeff)
{
    LOG("rescaling velocities with coefficient: " << coeff);
    foreach (particle& p, part) {
	p.v *= coeff;
    }
}

/**
 * initialize random number generator with seed
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::rng(unsigned int seed)
{
    rng_.set(seed);
    LOG("initializing random number generator with seed: " << seed);
}

/**
 * initialize random number generator from state
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::rng(gsl::gfsr4::state_type const& state)
{
    rng_.restore(state);
    LOG("restoring random number generator from state");
}

/**
 * place particles on a face-centered cubic (fcc) lattice
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::lattice()
{
    if (mixture_ == BINARY) {
	LOG("randomly placing A and B particles on fcc lattice");
	random_binary_types();
    }
    else {
	LOG("placing particles on fcc lattice");
    }

    // particles per 2- or 3-dimensional unit cell
    const unsigned int m = 2 * (dimension - 1);
    // lower boundary for number of particles per lattice dimension
    unsigned int n = static_cast<unsigned int>(std::pow(npart / m, 1. / dimension));
    // lower boundary for total number of lattice sites
    unsigned int N = m * static_cast<unsigned int>(pow(n, dimension));

    if (N < npart) {
	n += 1;
	N = m * static_cast<unsigned int>(pow(n, dimension));
    }
    if (N > npart) {
	LOG_WARNING("lattice not fully occupied (" << N << " sites)");
    }

    // lattice distance
    float_type a = box_ / n;
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

    // sort particles after binary mixture species for trajectory output
    struct compare
    {
	static bool _(particle const& p1, particle const& p2)
	{
	    return (p1.type < p2.type);
	}
    };
    std::stable_sort(this->part.begin(), this->part.end(), compare::_);

    // update cell lists
    update_cells();

    if (mixture_ == BINARY) {
	// update Verlet neighbour lists
	update_neighbours<true>();
	// calculate forces, potential energy and virial equation sum
	compute_forces<true>();
    }
    else {
	update_neighbours<false>();
	compute_forces<false>();
    }

    // reset sum over maximum velocity magnitudes to zero
    v_max_sum = 0;
}

/**
 * initialise velocities from Maxwell-Boltzmann distribution
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::temperature(double value)
{
    LOG("initializing velocities from Maxwell-Boltzmann distribution at temperature: " << value);

    // initialize force to zero for first leapfrog half step
    foreach (particle& p, part) {
	p.f = 0;
    }
    // initialize sum over maximum velocity magnitudes since last neighbour lists update
    v_max_sum = 0;

    boltzmann(value);
}

/**
 * set system temperature according to Maxwell-Boltzmann distribution
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::boltzmann(double temp)
{
    // center of mass velocity
    vector_type v_cm = 0;
    // mean squared velocity
    double vv = 0;

    // generate random Maxwell-Boltzmann distributed velocity
    foreach (particle& p, part) {
	rng_.gaussian(p.v, static_cast<float_type>(temp));
	v_cm += p.v;
    }
    v_cm /= npart;

    // set center of mass velocity to zero
    foreach (particle& p, part) {
	p.v -= v_cm;
	vv += p.v * p.v;
    }
    vv /= npart;

    // rescale velocities to accurate temperature
    double s = std::sqrt(temp * dimension / vv);
    foreach (particle& p, part) {
	p.v *= s;
    }
}

/**
 * randomly assign particles types in a binary mixture
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::random_binary_types()
{
    // create view on particle list
    std::vector<typename particle::ref> part;
    foreach (particle& p, this->part) {
	part.push_back(boost::ref(p));
    }

    // shuffle view and assign particles types
    rng_.shuffle(part);
    foreach (particle& p, range(part.begin(), part.begin() + mpart[0])) {
	p.type = particle::A;
    }
    foreach (particle& p, range(part.begin() + mpart[0], part.end())) {
	p.type = particle::B;
    }
}

/**
 * update cell lists
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::update_cells()
{
    // empty cell lists without memory reallocation
    foreach (cell_list& c, range(cell.data(), cell.data() + cell.num_elements())) {
	c.clear();
    }
    // add particles to cells
    foreach (particle& p, part) {
	compute_cell(p.r).push_back(boost::ref(p));
    }
}

/**
 * returns cell list which a particle belongs to
 */
template <int dimension>
typename ljfluid<ljfluid_impl_host, dimension>::cell_list&
ljfluid<ljfluid_impl_host, dimension>::compute_cell(vector_type r)
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
    cell_index index;
    for (int i = 0; i < dimension; ++i) {
	index[i] = (unsigned int)(r[i]) % ncell;
    }
    return cell(index);
}

/**
 * update neighbour lists
 */
template <int dimension>
template <bool binary>
void ljfluid<ljfluid_impl_host, dimension>::update_neighbours()
{
    cell_index i;
    for (i[0] = 0; i[0] < ncell; ++i[0]) {
	for (i[1] = 0; i[1] < ncell; ++i[1]) {
	    if (dimension == 3) {
		for (i[2] = 0; i[2] < ncell; ++i[2]) {
		    update_cell_neighbours<binary>(i);
		}
	    }
	    else {
		update_cell_neighbours<binary>(i);
	    }
	}
    }
}

/**
 * update neighbour lists for a single cell
 */
template <int dimension>
template <bool binary>
void ljfluid<ljfluid_impl_host, dimension>::update_cell_neighbours(cell_index const& i)
{
    foreach (particle& p, cell(i)) {
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
			compute_cell_neighbours<false, binary>(p, cell(k));
		    }
		}
		else {
		    // visit half of 8 neighbour cells due to pair potential
		    if (j[0] == 0 && j[1] == 0) {
			goto out;
		    }
		    // update neighbour list of particle
		    cell_index k;
		    for (int n = 0; n < dimension; ++n) {
			k[n] = (i[n] + ncell + j[n]) % ncell;
		    }
		    compute_cell_neighbours<false, binary>(p, cell(k));
		}
	    }
	}
out:
	// visit this cell
	compute_cell_neighbours<true, binary>(p, cell(i));
    }
}

/**
 * update neighbour list of particle
 */
template <int dimension>
template <bool same_cell, bool binary>
void ljfluid<ljfluid_impl_host, dimension>::compute_cell_neighbours(particle& p1, cell_list& c)
{
    foreach (particle& p2, c) {
	// skip identical particle and particle pair permutations if same cell
	if (same_cell && p2.tag <= p1.tag)
	    continue;

	// particle distance vector
	vector_type r = p1.r - p2.r;
	// binary particles type
	unsigned int type = (binary ? (p1.type + p2.type) : 0);
	// enforce periodic boundary conditions
	r -= round(r / box_) * box_;
	// squared particle distance
	float_type rr = r * r;

	// enforce cutoff radius with neighbour list skin
	if (rr >= static_cast<float_type>(rr_cut_skin[type]))
	    continue;

	// add particle to neighbour list
	p1.neighbour.push_back(boost::ref(p2));
    }
}

/**
 * compute Lennard-Jones forces
 */
template <int dimension>
template <bool binary>
void ljfluid<ljfluid_impl_host, dimension>::compute_forces()
{
    // initialize particle forces to zero
    foreach (particle& p, part) {
	p.f = 0;
    }

    // potential energy
    en_pot = 0;
    // virial equation sum
    virial.assign(binary ? 2 : 1, 0);

    foreach (particle& p1, part) {
	// calculate pairwise Lennard-Jones force with neighbour particles
	foreach (particle& p2, p1.neighbour) {
	    // particle distance vector
	    vector_type r = p1.r - p2.r;
	    // binary particles type
	    unsigned int type = (binary ? (p1.type + p2.type) : 0);
	    // enforce periodic boundary conditions
	    r -= round(r / box_) * box_;
	    // squared particle distance
	    float_type rr = r * r;

	    // enforce cutoff radius
	    if (rr >= rr_cut[type])
		continue;

	    // compute Lennard-Jones force in reduced units
	    float_type sigma2 = (binary ? sigma2_[type] : 1);
	    float_type eps = (binary ? epsilon_[type] : 1);
	    float_type rri = sigma2 / rr;
	    float_type r6i = rri * rri * rri;
	    float_type fval = 48 * rri * r6i * (r6i - 0.5) * (eps / sigma2);
	    float_type pot = (4 * r6i * (r6i - 1) - en_cut) * eps;

	    if (potential_ == C2POT) {
		compute_smooth_potential<binary>(std::sqrt(rr), fval, pot, type);
	    }

	    // add force contribution to both particles
	    p1.f += r * fval;
	    p2.f -= r * fval;

	    // add contribution to potential energy
	    en_pot += pot;

	    // add contribution to virial equation sum
	    float_type vir = 0.5 * rr * fval;
	    virial[p1.type][0] += vir;
	    virial[p2.type][0] += vir;

	    // compute off-diagonal virial stress tensor elements
	    if (dimension == 3) {
		vir = 0.5 * r[1] * r[2] * fval;
		virial[p1.type][1] += vir;
		virial[p2.type][1] += vir;

		vir = 0.5 * r[2] * r[0] * fval;
		virial[p1.type][2] += vir;
		virial[p2.type][2] += vir;

		vir = 0.5 * r[0] * r[1] * fval;
		virial[p1.type][3] += vir;
		virial[p2.type][3] += vir;
	    }
	    else {
		vir = 0.5 * r[0] * r[1] * fval;
		virial[p1.type][1] += vir;
		virial[p2.type][1] += vir;
	    }
	}
    }

    en_pot /= npart;

    // ensure that system is still in valid state
    if (std::isinf(en_pot)) {
	throw potential_energy_divergence();
    }
}

/**
 * compute C²-smooth potential
 */
template <int dimension>
template <bool binary>
void ljfluid<ljfluid_impl_host, dimension>::compute_smooth_potential(float_type r, float_type& fval, float_type& pot, unsigned int type)
{
    float_type y = r - r_cut[binary ? type : 0];
    float_type x2 = y * y * rri_smooth;
    float_type x4 = x2 * x2;
    float_type x4i = 1 / (1 + x4);
    // smoothing function
    float_type h0_r = x4 * x4i;
    // first derivative times (r_smooth)^(-1) [sic!]
    float_type h1_r = 4 * y * rri_smooth * x2 * x4i * x4i;
    // apply smoothing function to obtain C¹ force function
    fval = h0_r * fval - h1_r * (pot / r);
    // apply smoothing function to obtain C² potential function
    pot = h0_r * pot;
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::leapfrog_half()
{
    float_type vv_max = 0;

    foreach (particle& p, part) {
	// half step velocity
	p.v += p.f * (static_cast<float_type>(timestep_) / 2);
	// full step position
	p.r += p.v * static_cast<float_type>(timestep_);
	// maximum squared velocity
	vv_max = std::max(vv_max, p.v * p.v);
    }

    v_max_sum += std::sqrt(vv_max);
}

/**
 * second leapfrog step of integration of equations of motion
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::leapfrog_full()
{
    foreach (particle& p, part) {
	// full step velocity
	p.v += p.f * (static_cast<float_type>(timestep_) / 2);
    }
}

/**
 * MD simulation step
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::mdstep()
{
    // nanosecond resolution process times
    boost::array<timespec, 5> t;

    // calculate particle positions
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[0]);
    leapfrog_half();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[1]);

    if (v_max_sum * timestep_ > r_skin / 2) {
	// update cell lists
	update_cells();
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[2]);
	// update Verlet neighbour lists
	if (mixture_ == BINARY) {
	    update_neighbours<true>();
	}
	else {
	    update_neighbours<false>();
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[3]);
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0;

	m_times["update_cells"] += t[2] - t[1];
	m_times["update_neighbours"] += t[3] - t[2];
    }

    // calculate forces, potential energy and virial equation sum
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[2]);
    if (mixture_ == BINARY) {
	compute_forces<true>();
    }
    else {
	compute_forces<false>();
    }
    // calculate velocities
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[3]);
    if (thermostat_steps && ++thermostat_count > thermostat_steps) {
	boltzmann(thermostat_temp);
    }
    else {
	leapfrog_full();
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t[4]);

    if (thermostat_steps && thermostat_count > thermostat_steps) {
	// reset MD steps since last heatbath coupling
	thermostat_count = 0;
	m_times["boltzmann"] += t[4] - t[3];
	m_times["velocity_verlet"] += t[1] - t[0];
    }
    else {
	m_times["velocity_verlet"] += (t[1] - t[0]) + (t[4] - t[3]);
    }
    m_times["update_forces"] += t[3] - t[2];
    m_times["mdstep"] += t[4] - t[0];
}

template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::sample(host_sample_type& sample) const
{
    typedef typename host_sample_type::value_type sample_type;
    typedef typename sample_type::position_sample_vector position_sample_vector;
    typedef typename sample_type::position_sample_ptr position_sample_ptr;
    typedef typename sample_type::velocity_sample_vector velocity_sample_vector;
    typedef typename sample_type::velocity_sample_ptr velocity_sample_ptr;

    for (size_t n = 0, i = 0; n < npart; ++i) {
	// allocate memory for trajectory sample
	position_sample_ptr r(new position_sample_vector);
	velocity_sample_ptr v(new velocity_sample_vector);
	sample.push_back(sample_type(r, v));
	r->reserve(mpart[i]);
	v->reserve(mpart[i]);
	// assign particle positions and velocities of homogenous type
	for (size_t j = 0; j < mpart[i]; ++j, ++n) {
	    // periodically extended particle position
	    r->push_back(part[n].r);
	    // particle velocity
	    v->push_back(part[n].v);
	}
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::sample(energy_sample_type& sample) const
{
    typedef typename std::vector<particle>::const_iterator iterator;

    // virial tensor trace and off-diagonal elements for particle species
    sample.virial = virial;

    sample.vv = 0;
    sample.v_cm = 0;

    for (iterator p = part.begin(); p != part.end(); ++p) {
	// kinetic terms of virial stress tensor
	sample.virial[p->type][0] += p->v * p->v;
	if (dimension == 3) {
	    sample.virial[p->type][1] += p->v[1] * p->v[2];
	    sample.virial[p->type][2] += p->v[2] * p->v[0];
	    sample.virial[p->type][3] += p->v[0] * p->v[1];
	}
	else {
	    sample.virial[p->type][1] += p->v[0] * p->v[1];
	}

	sample.vv += p->v * p->v;
	sample.v_cm += p->v;
    }

    for (size_t i = 0; i < sample.virial.size(); ++i) {
	sample.virial[i] /= mpart[i];
    }

    // mean potential energy per particle
    sample.en_pot = en_pot;
    // mean squared velocity per particle
    sample.vv /= npart;
    // mean velocity per particle
    sample.v_cm /= npart;
}

/**
 * write parameters to HDF5 parameter group
 */
template <int dimension>
void ljfluid<ljfluid_impl_host, dimension>::param(H5param& param) const
{
    _Base::param(param);

    H5xx::group node(param["mdsim"]);
    node["cells"] = ncell;
    node["cell_length"] = cell_length_;
    node["neighbour_skin"] = r_skin;
}

} // namespace ljgpu

#undef foreach
#undef range

#endif /* ! LJGPU_MDSIM_LJFLUID_HOST_HPP */
