/* Lennard-Jones fluid simulation using CUDA
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
#include <cmath>
#include <cuda_wrapper.hpp>
#include <stdint.h>
#include "H5param.hpp"
#include "exception.hpp"
#include "log.hpp"
#include "statistics.hpp"
#include "rand48.hpp"
#include "vector2d.hpp"
#include "vector3d.hpp"
#include "gpu/ljfluid_glue.hpp"


#define foreach BOOST_FOREACH

namespace mdsim
{

/**
 * Lennard-Jones fluid simulation using CUDA
 */
template <unsigned dimension, typename T>
class ljfluid
{
public:
    /** device floating-point vector type */
    typedef typename cuda::types::vector<dimension, typename T::value_type>::type floatn;

public:
    /** initialize fixed simulation parameters */
    ljfluid();
    /** set number of particles in system */
    void particles(unsigned int value);
    /** set particle density */
    void density(float value);
    /** set periodic box length */
    void box(float value);
#ifdef USE_CELL
    /** set desired average cell occupancy */
    void cell_occupancy(float value);
#endif
    /** set number of CUDA execution threads */
    void threads(unsigned int value);
    /** restore system state from phase space sample */
    template <typename V> void restore(V visitor);

    /** seed random number generator */
    void rng(unsigned int seed);
    /** restore random number generator from state */
    void rng(mdsim::rand48::state_type const& state);
    /** place particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(float temp);
    /** set simulation timestep */
    void timestep(float value);

    /** get number of particles */
    unsigned int const& particles() const;
#ifdef USE_CELL
    /** get number of cells per dimension */
    unsigned int const& cells() const;
    /** get total number of cell placeholders */
    unsigned int const& placeholders() const;
    /** get cell length */
    float const& cell_length() const;
    /** get effective average cell occupancy */
    float const& cell_occupancy() const;
    /** get number of placeholders per cell */
    unsigned int const& cell_size() const;
#endif
    /** get number of CUDA execution blocks */
    unsigned int blocks() const;
    /** get number of CUDA execution threads */
    unsigned int threads() const;
    /** get particle density */
    float const& density() const;
    /** get periodic box length */
    float const& box() const;
    /** get simulation timestep */
    float const& timestep() const;
    /** get potential cutoff distance */
    float const& cutoff_distance() const;

    /** stream MD simulation step on GPU */
    void mdstep();
    /** copy MD simulation step results from GPU to host */
    void sample();
    /** sample trajectory */
    template <typename V> void sample(V visitor) const;

private:
    /** number of particles in system */
    unsigned int npart;
#ifdef USE_CELL
    /** number of cells per dimension */
    unsigned int ncell;
    /** total number of cell placeholders */
    unsigned int nplace;
    /** cell length */
    float cell_length_;
    /** effective average cell occupancy */
    float cell_occupancy_;
    /** number of placeholders per cell */
    unsigned int cell_size_;
#endif
    /** particle density */
    float density_;
    /** periodic box length */
    float box_;
    /** simulation timestep */
    float timestep_;
    /** cutoff distance for shifted Lennard-Jones potential */
    float r_cut;

    /** system state in page-locked host memory */
    struct {
	/** periodically reduced particle positions */
	cuda::host::vector<T> r;
	/** periodically extended particle positions */
	cuda::host::vector<T> R;
	/** particle velocities */
	cuda::host::vector<T> v;
#ifndef USE_CELL
	/** potential energies per particle */
	cuda::host::vector<float> en;
	/** virial equation sums per particle */
	cuda::host::vector<float> virial;
#endif
    } h_part;

    /** mean potential energy per particle */
    float en_pot_;
    /** mean virial equation sum per particle */
    float virial_;

#ifdef USE_CELL
    /** cell placeholders in page-locked host memory */
    struct {
	/** periodically reduced particle positions */
	cuda::host::vector<T> r;
	/** periodically extended particle positions */
	cuda::host::vector<T> R;
	/** particle velocities */
	cuda::host::vector<T> v;
	/** particle number tags */
	cuda::host::vector<int> n;
	/** potential energies per particle */
	cuda::host::vector<float> en;
	/** virial equation sums per particle */
	cuda::host::vector<float> virial;
    } h_cell;
#endif

#ifdef USE_CELL
    /** system state in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<floatn> r;
	/** periodically extended particle positions */
	cuda::vector<floatn> R;
	/** particle velocities */
	cuda::vector<floatn> v;
	/** particle number tags */
	cuda::vector<int> n;
	/** particle forces */
	cuda::vector<floatn> f;
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
    } g_cell;

    /** system state double buffer in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<floatn> r;
	/** periodically extended particle positions */
	cuda::vector<floatn> R;
	/** particle velocities */
	cuda::vector<floatn> v;
	/** particle number tags */
	cuda::vector<int> n;
    } g_cell2;
#else
    /** system state in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<floatn> r;
	/** periodically extended particle positions */
	cuda::vector<floatn> R;
	/** particle velocities */
	cuda::vector<floatn> v;
	/** particle forces */
	cuda::vector<floatn> f;
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
    } g_part;
#endif

    /** random number generator */
    mdsim::rand48 rng_;
    /** CUDA execution dimensions */
    cuda::config dim_;
#ifdef USE_CELL
    /** CUDA execution dimensions for cell-specific kernels */
    cuda::config dim_cell_;
#endif
    /** CUDA asynchronous execution */
    cuda::stream stream_;
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
    float rr_cut = r_cut * r_cut;
    // potential energy at cutoff distance
    float rri_cut = 1. / rr_cut;
    float r6i_cut = rri_cut * rri_cut * rri_cut;
    float en_cut = 2. * r6i_cut * (r6i_cut - 1.);

    try {
	cuda::copy(rr_cut, gpu::ljfluid::rr_cut);
	cuda::copy(en_cut, gpu::ljfluid::en_cut);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy cutoff parameters to device symbols");
    }
}

/**
 * set number of particles in system
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::particles(unsigned int value)
{
    // validate particle number
    if (value < 1) {
	throw exception("invalid number of particles");
    }
    // set particle number
    npart = value;
    LOG("number of particles: " << npart);
    // copy particle number to device symbol
    try {
	cuda::copy(npart, gpu::ljfluid::npart);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy particle number to device symbol");
    }

#ifndef USE_CELL
    // allocate global device memory for system state
    try {
	g_part.r.resize(npart);
	g_part.R.resize(npart);
	g_part.v.resize(npart);
	g_part.f.resize(npart);
	g_part.en.resize(npart);
	g_part.virial.resize(npart);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for system state");
    }
#endif

    // allocate page-locked host memory for system state
    try {
	h_part.r.resize(npart);
	h_part.R.resize(npart);
	h_part.v.resize(npart);
	// particle forces reside only in GPU memory
#ifndef USE_CELL
	h_part.en.resize(npart);
	h_part.virial.resize(npart);
#endif
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate page-locked host memory for system state");
    }

#ifndef USE_CELL
    // set particle forces to zero for first integration of differential equations of motion
    try {
	fill(h_part.v.begin(), h_part.v.end(), 0.);
	cuda::copy(h_part.v, g_part.f, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to initialize forces to zero");
    }
#endif
}

/**
 * set particle density
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::density(float value)
{
    // set particle density
    density_ = value;
    LOG("particle density: " << density_);

    // compute periodic box length
    box_ = powf(npart / density_, 1. / dimension);
    LOG("periodic simulation box length: " << box_);
    // copy periodic box length to device symbol
    try {
	cuda::copy(box_, gpu::ljfluid::box);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy periodic box length to device symbol");
    }
}

/**
 * set periodic box length
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::box(float value)
{
    // set periodic box length
    box_ = value;
    LOG("periodic simulation box length: " << box_);
    // copy periodic box length to device symbol
    try {
	cuda::copy(box_, gpu::ljfluid::box);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy periodic box length to device symbol");
    }

    // compute particle density
    density_ = npart / powf(box_, dimension);
    LOG("particle density: " << density_);
}

#ifdef USE_CELL
/**
 * set desired average cell occupancy
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::cell_occupancy(float value)
{
    LOG("desired average cell occupancy: " << value);

    // fixed cell size due to fixed number of CUDA execution threads per block
    cell_size_ = CELL_SIZE;

    // optimal cell length to yield desired average cell occupancy
    cell_length_ = std::pow((cell_size_ * value) / density_, 1.f / dimension);

    // set number of cells per dimension
    ncell = std::floor(box_ / std::max(r_cut, cell_length_));
    LOG("number of cells per dimension: " << ncell);

    if (ncell < 3) {
	throw exception("number of cells per dimension must be at least 3");
    }

    // set cell length from integer number of cells
    cell_length_ = box_ / ncell;
    LOG("cell length: " << cell_length_);

    // set total number of cell placeholders
    nplace = pow(ncell, dimension) * cell_size_;
    LOG("number of cell placeholders: " << nplace);

    // set effective average cell occupancy
    cell_occupancy_ = npart * 1. / nplace;
    LOG("effective average cell occupancy: " << cell_occupancy_);

    if (cell_occupancy_ > 1.) {
	throw exception("average cell occupancy must not be larger than 1.0");
    }
    else if (cell_occupancy_ > 0.5) {
	LOG_WARNING("average cell occupancy is larger than 0.5");
    }

    // copy cell parameters to device symbols
    try {
	cuda::copy(ncell, gpu::ljfluid::ncell);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy cell parameters to device symbols");
    }
}
#endif /* USE_CELL */

/**
 * set number of CUDA execution threads
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::threads(unsigned int value)
{
    // query CUDA device properties
    cuda::device::properties prop;
    try {
	prop = cuda::device::properties(cuda::device::get());
    }
    catch (cuda::error const& e) {
	throw exception("failed to query CUDA device properties");
    }

    // validate number of CUDA execution threads
    if (value < 1) {
	throw exception("invalid number of CUDA execution threads");
    }
    if (value > prop.max_threads_per_block()) {
	throw exception("number of CUDA execution threads exceeds maximum number of threads per block");
    }
    if (value & (value - 1)) {
	LOG_WARNING("number of CUDA execution threads not a power of 2");
    }
    if (value % prop.warp_size()) {
	LOG_WARNING("number of CUDA execution threads not a multiple of warp size (" << prop.warp_size() << ")");
    }

    // set CUDA execution dimensions
    dim_ = cuda::config(dim3((npart + value - 1) / value), dim3(value));
    LOG("number of CUDA execution blocks: " << dim_.blocks_per_grid());
    LOG("number of CUDA execution threads: " << dim_.threads_per_block());

    if (dim_.threads() != npart) {
	LOG_WARNING("number of particles (" << npart << ") not a multiple of number of CUDA execution threads (" << dim_.threads() << ")");
    }

#ifdef USE_CELL
    // set CUDA execution dimensions for cell-specific kernels
    dim_cell_ = cuda::config(dim3(powf(ncell, dimension)), dim3(cell_size_));
    LOG("number of cell CUDA execution blocks: " << dim_cell_.blocks_per_grid());
    LOG("number of cell CUDA execution threads: " << dim_cell_.threads_per_block());

    // allocate page-locked host memory for cell placeholders
    try {
	h_cell.r.resize(dim_cell_.threads());
	h_cell.R.resize(dim_cell_.threads());
	h_cell.v.resize(dim_cell_.threads());
	h_cell.n.resize(dim_cell_.threads());
	// particle forces reside only in GPU memory
	h_cell.en.resize(dim_cell_.threads());
	h_cell.virial.resize(dim_cell_.threads());
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate page-locked host memory cell placeholders");
    }

    // allocate global device memory for cell placeholders
    try {
	g_cell.r.resize(dim_cell_.threads());
	g_cell.R.resize(dim_cell_.threads());
	g_cell.v.resize(dim_cell_.threads());
	g_cell.n.resize(dim_cell_.threads());
	g_cell.f.resize(dim_cell_.threads());
	g_cell.en.resize(dim_cell_.threads());
	g_cell.virial.resize(dim_cell_.threads());
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory cell placeholders");
    }

    // allocate global device memory for double buffers
    try {
	g_cell2.r.resize(dim_cell_.threads());
	g_cell2.R.resize(dim_cell_.threads());
	g_cell2.v.resize(dim_cell_.threads());
	g_cell2.n.resize(dim_cell_.threads());
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory double buffers");
    }

    // set particle forces to zero for first integration of differential equations of motion
    try {
	fill(h_cell.v.begin(), h_cell.v.end(), 0.);
	cuda::copy(h_cell.v, g_cell.f, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to initialize forces to zero");
    }
#else /* USE_CELL */
    // allocate global device memory for placeholder particles
    try {
	g_part.r.reserve(dim_.threads());
	g_part.R.reserve(dim_.threads());
	g_part.v.reserve(dim_.threads());
	g_part.f.reserve(dim_.threads());
	g_part.en.reserve(dim_.threads());
	g_part.virial.reserve(dim_.threads());
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for placeholder particles");
    }
#endif /* USE_CELL */

    // change random number generator dimensions
    try {
	rng_.resize(dim_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to change random number generator dimensions");
    }
}

/**
 * restore system state from phase space sample
 */
template <unsigned dimension, typename T>
template <typename V>
void ljfluid<dimension, T>::restore(V visitor)
{
    // read phase space sample
    visitor(h_part.r, h_part.v);

    try {
#ifdef USE_CELL
	// copy periodically reduced particle positions from host to GPU
	cuda::vector<floatn> g_r(npart);
	cuda::copy(h_part.r, g_r, stream_);

	// assign particles to cells
	gpu::ljfluid::assign_cells.configure(dim_cell_, stream_);
	gpu::ljfluid::assign_cells(g_r.data(), g_cell.r.data(), g_cell.n.data());
	// copy particle number tags from GPU to host
	cuda::copy(g_cell.n, h_cell.n, stream_);
	// copy particle positions from GPU to host
	cuda::copy(g_cell.r, h_cell.r, stream_);
	// replicate to periodically extended positions
	cuda::copy(g_cell.r, g_cell.R, stream_);
	cuda::copy(h_cell.r, h_cell.R, stream_);

	// recalculate forces for first leapfrog half step
	gpu::ljfluid::mdstep.configure(dim_cell_, stream_);
	gpu::ljfluid::mdstep(g_cell.r.data(), g_cell.v.data(), g_cell.f.data(), g_cell.n.data(), g_cell.en.data(), g_cell.virial.data());

	// assign velocities to cell placeholders
	for (unsigned int i = 0; i < nplace; ++i) {
	    if (IS_REAL_PARTICLE(h_cell.n[i])) {
		h_cell.v[i] = h_part.v[h_cell.n[i]];
	    }
	}
	// copy particle velocities from host to GPU (after force calculation!)
	cuda::copy(h_cell.v, g_cell.v, stream_);
#else
	// copy periodically reduced particle positions from host to GPU
	cuda::copy(h_part.r, g_part.r, stream_);
	// replicate to periodically extended particle positions
	cuda::copy(g_part.r, g_part.R, stream_);
	cuda::copy(h_part.r, h_part.R, stream_);
	// recalculate forces for first leapfrog half step
	gpu::ljfluid::mdstep.configure(dim_, dim_.threads_per_block() * sizeof(T), stream_);
	gpu::ljfluid::mdstep(g_part.r.data(), g_part.v.data(), g_part.f.data(), g_part.en.data(), g_part.virial.data());
	// copy particle velocities from host to GPU (after force calculation!)
	cuda::copy(h_part.v, g_part.v, stream_);
#endif

	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to restore system state from phase space sample");
    }
}

/**
 * seed random number generator
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::rng(unsigned int seed)
{
    LOG("random number generator seed: " << seed);
    try {
	rng_.set(seed);
    }
    catch (cuda::error const& e) {
	throw exception("failed to seed random number generator");
    }
}

/**
 * restore random number generator from state
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::rng(mdsim::rand48::state_type const& state)
{
    try {
	rng_.restore(state);
    }
    catch (cuda::error const& e) {
	throw exception("failed to restore random number generator state");
    }
}

/**
 * place particles on a face-centered cubic (fcc) lattice
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::lattice()
{
    LOG("placing particles on face-centered cubic (fcc) lattice");
    try {
#ifdef USE_CELL
	cuda::vector<floatn> g_r(npart);
	g_r.reserve(dim_.threads());
	// compute particle lattice positions on GPU
	gpu::ljfluid::lattice.configure(dim_, stream_);
	gpu::ljfluid::lattice(g_r.data());
	// copy particle positions from GPU to host
	cuda::copy(g_r, h_part.r, stream_);
	// assign particles to cells
	gpu::ljfluid::assign_cells.configure(dim_cell_, stream_);
	gpu::ljfluid::assign_cells(g_r.data(), g_cell.r.data(), g_cell.n.data());
	// copy particle number tags from GPU to host
	cuda::copy(g_cell.n, h_cell.n, stream_);
	// copy particle positions from GPU to host
	cuda::copy(g_cell.r, h_cell.r, stream_);
	// copy particle positions to periodically extended positions
	cuda::copy(g_cell.r, g_cell.R, stream_);
	cuda::copy(h_cell.r, h_cell.R, stream_);
#else
	// compute particle lattice positions on GPU
	gpu::ljfluid::lattice.configure(dim_, stream_);
	gpu::ljfluid::lattice(g_part.r.data());
	// copy particle positions from GPU to host
	cuda::copy(g_part.r, h_part.r, stream_);
	// copy particle positions to periodically extended positions
	cuda::copy(g_part.r, g_part.R, stream_);
	cuda::copy(h_part.r, h_part.R, stream_);
#endif

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute particle lattice positions on GPU");
    }
}

/**
 * set system temperature according to Maxwell-Boltzmann distribution
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::temperature(float temp)
{
    LOG("initializing velocities from Maxwell-Boltzmann distribution at temperature: " << temp);
    try {
#ifdef USE_CELL
	cuda::vector<floatn> g_v(npart);
	g_v.reserve(dim_.threads());
	// set velocities using Maxwell-Boltzmann distribution at temperature
	gpu::ljfluid::boltzmann.configure(dim_, stream_);
	gpu::ljfluid::boltzmann(g_v.data(), temp, rng_.data());
	// copy particle velocities from GPU to host
	cuda::copy(g_v, h_part.v, stream_);
#else
	// set velocities using Maxwell-Boltzmann distribution at temperature
	gpu::ljfluid::boltzmann.configure(dim_, stream_);
	gpu::ljfluid::boltzmann(g_part.v.data(), temp, rng_.data());
	// copy particle velocities from GPU to host
	cuda::copy(g_part.v, h_part.v, stream_);
#endif

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute Maxwell-Boltzmann distributed velocities on GPU");
    }

    // compute center of mass velocity
    T v_cm = mean(h_part.v.begin(), h_part.v.end());
    // set center of mass velocity to zero
    foreach (T& v, h_part.v) {
	v -= v_cm;
    }

    try {
#ifdef USE_CELL
	// assign velocities to cell placeholders
	for (unsigned int i = 0; i < nplace; ++i) {
	    if (IS_REAL_PARTICLE(h_cell.n[i])) {
		h_cell.v[i] = h_part.v[h_cell.n[i]];
	    }
	}
	// copy particle velocities from host to GPU
	cuda::copy(h_cell.v, g_cell.v, stream_);
#else
	// copy particle velocities from host to GPU
	cuda::copy(h_part.v, g_part.v, stream_);
#endif
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to set center of mass velocity to zero");
    }
}

/**
 * set simulation timestep
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::timestep(float value)
{
    // set simulation timestep
    timestep_ = value;
    LOG("simulation timestep: " << timestep_);
    // copy simulation timestep to device symbol
    try {
	cuda::copy(timestep_, gpu::ljfluid::timestep);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy simulation timestep to device symbol");
    }
}

/**
 * get number of particles
 */
template <unsigned dimension, typename T>
unsigned int const& ljfluid<dimension, T>::particles() const
{
    return npart;
}

#ifdef USE_CELL
/**
 * get number of cells per dimension
 */
template <unsigned dimension, typename T>
unsigned int const& ljfluid<dimension, T>::cells() const
{
    return ncell;
}

/**
 * get total number of cell placeholders
 */
template <unsigned dimension, typename T>
unsigned int const& ljfluid<dimension, T>::placeholders() const
{
    return nplace;
}

/**
 * get cell length
 */
template <unsigned dimension, typename T>
float const& ljfluid<dimension, T>::cell_length() const
{
    return cell_length_;
}

/**
 * get effective average cell occupancy
 */
template <unsigned dimension, typename T>
float const& ljfluid<dimension, T>::cell_occupancy() const
{
    return cell_occupancy_;
}

/**
 * get number of placeholders per cell
 */
template <unsigned dimension, typename T>
unsigned int const& ljfluid<dimension, T>::cell_size() const
{
    return cell_size_;
}
#endif

/**
 * get number of CUDA execution blocks
 */
template <unsigned dimension, typename T>
unsigned int ljfluid<dimension, T>::blocks() const
{
    return dim_.blocks_per_grid();
}

/**
 * get number of particles
 */
template <unsigned dimension, typename T>
unsigned int ljfluid<dimension, T>::threads() const
{
    return dim_.threads_per_block();
}

/**
 * get particle density
 */
template <unsigned dimension, typename T>
float const& ljfluid<dimension, T>::density() const
{
    return density_;
}

/**
 * get periodic box length
 */
template <unsigned dimension, typename T>
float const& ljfluid<dimension, T>::box() const
{
    return box_;
}

/**
 * get simulation timestep
 */
template <unsigned dimension, typename T>
float const& ljfluid<dimension, T>::timestep() const
{
    return timestep_;
}

/**
 * get potential cutoff distance
 */
template <unsigned dimension, typename T>
float const& ljfluid<dimension, T>::cutoff_distance() const
{
    return r_cut;
}

/**
 * stream MD simulation step on GPU
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::mdstep()
{
#ifdef USE_CELL
    // first leapfrog step of integration of differential equations of motion
    try {
	gpu::ljfluid::inteq.configure(dim_cell_, stream_);
	gpu::ljfluid::inteq(g_cell.r.data(), g_cell.R.data(), g_cell.v.data(), g_cell.f.data());
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream first leapfrog step on GPU");
    }

    // update cell lists
    try {
	gpu::ljfluid::update_cells.configure(dim_cell_, stream_);
	gpu::ljfluid::update_cells(g_cell.r.data(), g_cell.R.data(), g_cell.v.data(), g_cell.n.data(), g_cell2.r.data(), g_cell2.R.data(), g_cell2.v.data(), g_cell2.n.data());

	cuda::copy(g_cell2.r, g_cell.r, stream_);
	cuda::copy(g_cell2.R, g_cell.R, stream_);
	cuda::copy(g_cell2.v, g_cell.v, stream_);
	cuda::copy(g_cell2.n, g_cell.n, stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream cell list update on GPU");
    }

    // Lennard-Jones force calculation
    try {
	gpu::ljfluid::mdstep.configure(dim_cell_, stream_);
	gpu::ljfluid::mdstep(g_cell.r.data(), g_cell.v.data(), g_cell.f.data(), g_cell.n.data(), g_cell.en.data(), g_cell.virial.data());
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream force calculation on GPU");
    }
#else
    // first leapfrog step of integration of differential equations of motion
    try {
	gpu::ljfluid::inteq.configure(dim_, stream_);
	gpu::ljfluid::inteq(g_part.r.data(), g_part.R.data(), g_part.v.data(), g_part.f.data());
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream first leapfrog step on GPU");
    }

    // Lennard-Jones force calculation
    try {
	gpu::ljfluid::mdstep.configure(dim_, dim_.threads_per_block() * sizeof(T), stream_);
	gpu::ljfluid::mdstep(g_part.r.data(), g_part.v.data(), g_part.f.data(), g_part.en.data(), g_part.virial.data());
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream force calculation on GPU");
    }
#endif
}

/**
 * copy MD simulation step results from GPU to host
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::sample()
{
    // wait for MD simulation step on GPU to finish
    try {
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("MD simulation step on GPU failed");
    }

#ifdef USE_CELL
    // copy MD simulation step results from GPU to host
    try {
	// copy periodically reduce particles positions
	cuda::copy(g_cell.r, h_cell.r, stream_);
	// copy periodically extended particles positions
	cuda::copy(g_cell.R, h_cell.R, stream_);
	// copy particle velocities
	cuda::copy(g_cell.v, h_cell.v, stream_);
	// copy particle number tags
	cuda::copy(g_cell.n, h_cell.n, stream_);
	// copy potential energies per particle and virial equation sums
	cuda::copy(g_cell.en, h_cell.en, stream_);
	// copy virial equation sums per particle
	cuda::copy(g_cell.virial, h_cell.virial, stream_);

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy MD simulation step results from GPU to host");
    }

    // mean potential energy per particle
    en_pot_ = 0.;
    // mean virial equation sum per particle
    virial_ = 0.;
    // number of particles found in cells
    unsigned int count = 0;

    for (unsigned int i = 0; i < nplace; ++i) {
	// check if real particle
	if (IS_REAL_PARTICLE(h_cell.n[i])) {
	    // copy periodically reduced particle positions
	    h_part.r[h_cell.n[i]] = h_cell.r[i];
	    // copy periodically extended particle positions
	    h_part.R[h_cell.n[i]] = h_cell.R[i];
	    // copy particle velocities
	    h_part.v[h_cell.n[i]] = h_cell.v[i];

	    count++;

	    // calculate mean potential energy per particle
	    en_pot_ += (h_cell.en[i] - en_pot_) / count;
	    // calculate mean virial equation sum per particle
	    virial_ += (h_cell.virial[i] - virial_) / count;
	}
    }
    // validate number of particles
    if (count != npart) {
	throw exception("particle loss while updating cell lists");
    }
#else
    // copy MD simulation step results from GPU to host
    try {
	// copy periodically reduce particles positions
	cuda::copy(g_part.r, h_part.r, stream_);
	// copy periodically extended particles positions
	cuda::copy(g_part.R, h_part.R, stream_);
	// copy particle velocities
	cuda::copy(g_part.v, h_part.v, stream_);
	// copy potential energies per particle and virial equation sums
	cuda::copy(g_part.en, h_part.en, stream_);
	// copy virial equation sums per particle
	cuda::copy(g_part.virial, h_part.virial, stream_);

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy MD simulation step results from GPU to host");
    }

    // calculate mean potential energy per particle
    en_pot_ = mean(h_part.en.begin(), h_part.en.end());
    // calculate mean virial equation sum per particle
    virial_ = mean(h_part.virial.begin(), h_part.virial.end());
#endif

    // ensure that system is still in valid state after MD step
    if (std::isnan(en_pot_)) {
	throw exception("potential energy diverged due to excessive timestep or density");
    }
}

/**
 * sample trajectory
 */
template <unsigned dimension, typename T>
template <typename V>
void ljfluid<dimension, T>::sample(V visitor) const
{
    visitor(h_part.r, h_part.R, h_part.v, en_pot_, virial_);
}

} // namespace mdsim

#undef foreach

#endif /* ! MDSIM_LJFLUID_HPP */
