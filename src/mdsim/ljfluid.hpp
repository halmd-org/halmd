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
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <cmath>
#include <cuda_wrapper.hpp>
#include <stdint.h>
#include "H5param.hpp"
#include "H5xx.hpp"
#include "accumulator.hpp"
#include "exception.hpp"
#include "gpu/ljfluid_glue.hpp"
#include "log.hpp"
#include "perf.hpp"
#include "radix.hpp"
#include "rand48.hpp"
#include "sample.hpp"
#include "statistics.hpp"

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
#ifdef DIM_3D
    /** coalesced GPU single-precision floating-point types */
    typedef float4 g_vector;
#else
    typedef float2 g_vector;
#endif

public:
    /** initialise fixed simulation parameters */
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
    unsigned int const& particles() const { return npart; }
#ifdef USE_CELL
    /** get number of cells per dimension */
    unsigned int const& cells() const { return ncell; }
    /** get total number of cell placeholders */
    unsigned int const& placeholders() const { return nplace; }
    /** get cell length */
    float const& cell_length() const { return cell_length_; }
    /** get effective average cell occupancy */
    float const& cell_occupancy() const { return cell_occupancy_; }
    /** get number of placeholders per cell */
    unsigned int const& cell_size() const { return cell_size_; }
#endif
    /** get number of CUDA execution blocks */
    unsigned int blocks() const { return dim_.blocks_per_grid(); }
    /** get number of CUDA execution threads */
    unsigned int threads() const { return dim_.threads_per_block(); }
    /** get particle density */
    float const& density() const { return density_; }
    /** get periodic box length */
    float const& box() const { return box_; }
    /** get simulation timestep */
    float const& timestep() const { return timestep_; }
    /** get potential cutoff distance */
    float const& cutoff_distance() const { return r_cut; }
    /** returns and resets CUDA time statistics */
    perf_counters times();

    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;

    /** stream MD simulation step on GPU */
    void mdstep();
    /** synchronize MD simulation step on GPU */
    void synchronize();
    /** copy MD simulation step results from GPU to host */
    void sample();
    /** get trajectory sample */
    trajectory_sample<T> const& trajectory() const { return h_sample; }

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
    /** maximum velocity magnitude after last MD step */
    float v_max;
#ifdef USE_CELL
    /** cell skin */
    float r_skin;
    /** sum over maximum velocity magnitudes since last cell lists update */
    float v_max_sum;
#endif
#ifdef USE_SMOOTH_POTENTIAL
    /** potential smoothing function scale parameter */
    float r_smooth;
#endif

    /** trajectory sample in swappable host memory */
    trajectory_sample<T> h_sample;

#ifdef USE_CELL
    /** cell placeholders in page-locked host memory */
    struct {
	/** periodically reduced particle positions */
	cuda::host::vector<g_vector> r;
	/** periodically extended particle positions */
	cuda::host::vector<g_vector> R;
	/** particle velocities */
	cuda::host::vector<g_vector> v;
	/** particle number tags */
	cuda::host::vector<int> n;
	/** potential energies per particle */
	cuda::host::vector<float> en;
	/** virial equation sums per particle */
	cuda::host::vector<float> virial;
    } h_cell;
#else
    /** system state in page-locked host memory */
    struct {
	/** periodically reduced particle positions */
	cuda::host::vector<g_vector> r;
	/** periodically extended particle positions */
	cuda::host::vector<g_vector> R;
	/** particle velocities */
	cuda::host::vector<g_vector> v;
	/** potential energies per particle */
	cuda::host::vector<float> en;
	/** virial equation sums per particle */
	cuda::host::vector<float> virial;
    } h_part;
#endif

#ifdef USE_CELL
    /** system state in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<g_vector> r;
	/** periodically extended particle positions */
	cuda::vector<g_vector> R;
	/** particle velocities */
	cuda::vector<g_vector> v;
	/** particle number tags */
	cuda::vector<int> n;
	/** particle forces */
	cuda::vector<g_vector> f;
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
    } g_cell;

    /** system state double buffer in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<g_vector> r;
	/** periodically extended particle positions */
	cuda::vector<g_vector> R;
	/** particle velocities */
	cuda::vector<g_vector> v;
	/** particle number tags */
	cuda::vector<int> n;
    } g_cell2;
#else
    /** system state in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<g_vector> r;
	/** periodically extended particle positions */
	cuda::vector<g_vector> R;
	/** particle velocities */
	cuda::vector<g_vector> v;
	/** particle forces */
	cuda::vector<g_vector> f;
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
    } g_part;
#endif

    /** GPU random number generator */
    mdsim::rand48 rng_;
    /** GPU mdsim::radix sort for particle positions */
    mdsim::radix_sort<g_vector> radix_sort;
    /** CUDA execution dimensions */
    cuda::config dim_;
#ifdef USE_CELL
    /** CUDA execution dimensions for cell-specific kernels */
    cuda::config dim_cell_;
#endif
    /** CUDA asynchronous execution */
    cuda::stream stream_;
    /** CUDA events for kernel timing */
#ifdef USE_CELL
    boost::array<cuda::event, 5> event_;
#else
    boost::array<cuda::event, 3> event_;
#endif
    /** CUDA time statistics */
    perf_counters m_times;
};


/**
 * initialise fixed simulation parameters
 */
template <unsigned dimension, typename T>
ljfluid<dimension, T>::ljfluid()
{
    // suppress attractive tail of Lennard-Jones potential
    r_cut = std::pow(2, 1 / 6.f);
    LOG("potential cutoff distance: " << r_cut);

    // squared cutoff distance
    float rr_cut = std::pow(r_cut, 2);
    // potential energy at cutoff distance
    float rri_cut = 1 / rr_cut;
    float r6i_cut = rri_cut * rri_cut * rri_cut;
    float en_cut = 4 * r6i_cut * (r6i_cut - 1);

    LOG("potential cutoff energy: " << en_cut);

    try {
	cuda::copy(r_cut, gpu::ljfluid::r_cut);
	cuda::copy(rr_cut, gpu::ljfluid::rr_cut);
	cuda::copy(en_cut, gpu::ljfluid::en_cut);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy potential cutoff parameters to device symbols");
    }

#ifdef USE_SMOOTH_POTENTIAL
    // compute potential smoothing function scale parameter
    r_smooth = 0.001;
    LOG("potential smoothing function scale parameter: " << r_smooth);

    // squared inverse potential smoothing function scale parameter
    float rri_smooth = std::pow(r_smooth, -2);

    try {
	cuda::copy(rri_smooth, gpu::ljfluid::rri_smooth);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy potential smoothing function scale parameter to device symbol");
    }
#endif
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

    // allocate swappable host memory for trajectory sample
    try {
	h_sample.r.resize(npart);
	h_sample.R.resize(npart);
	h_sample.v.resize(npart);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate swappable host memory for trajectory sample");
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

    // allocate page-locked host memory for system state
    try {
	h_part.r.resize(npart);
	h_part.R.resize(npart);
	h_part.v.resize(npart);
	// particle forces reside only in GPU memory
	h_part.en.resize(npart);
	h_part.virial.resize(npart);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate page-locked host memory for system state");
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

    // optimal number of cells with given cell occupancy as upper boundary
    ncell = std::ceil(std::pow(npart / (value * cell_size_), 1.f / dimension));

    // set number of cells per dimension, respecting cutoff distance
    ncell = std::min(ncell, uint(box_ / r_cut));
    LOG("number of cells per dimension: " << ncell);

    if (ncell < 3) {
	throw exception("number of cells per dimension must be at least 3");
    }

    // derive cell length from number of cells
    cell_length_ = box_ / ncell;
    LOG("cell length: " << cell_length_);
    // set cell skin
    r_skin = std::max(0.f, cell_length_ - r_cut);
    LOG("cell skin: " << r_skin);

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

#ifndef USE_CELL
    if (dim_.threads() != npart) {
	LOG_WARNING("number of particles (" << npart << ") not a multiple of number of CUDA execution threads (" << dim_.threads() << ")");
    }
#endif

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

    try {
	// compute optimal number of blocks for GeForce 8800 with 16 multiprocessors
	const uint threads = dim_.threads_per_block();
	const uint max_blocks = (16 * 512) / (threads * gpu::radix::BUCKETS_PER_THREAD / 2);
	const uint blocks = std::min((npart + 2 * threads - 1) / (2 * threads), max_blocks);

	LOG("number of CUDA blocks for radix sort: " << blocks);
	LOG("number of CUDA threads for radix sort: " << threads);

	// allocate global device memory for radix sort
	radix_sort.resize(npart, blocks, threads);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for radix sort");
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
    visitor(h_sample.r, h_sample.v);

    try {
#ifdef USE_CELL
	// copy periodically reduced particle positions from host to GPU
	for (unsigned int i = 0; i < npart; ++i) {
	    h_cell.r[i] = make_float(h_sample.r[i]);
	}
	cuda::copy(h_cell.r, g_cell2.r, stream_);
	// assign particles to cells
	event_[0].record(stream_);
	cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	gpu::ljfluid::assign_cells(g_cell2.r, g_cell.r, g_cell.n);
	event_[1].record(stream_);
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0;
	// replicate particle positions to periodically extended positions
	cuda::copy(g_cell.r, g_cell.R, stream_);
	// calculate forces, potential energy and virial equation sum
	cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	gpu::ljfluid::mdstep(g_cell.r, g_cell.v, g_cell.f, g_cell.n, g_cell.en, g_cell.virial);

	// copy particle number tags from GPU to host
	cuda::copy(g_cell.n, h_cell.n, stream_);
	// wait for CUDA operations to finish
	stream_.synchronize();

	// maximum squared velocity
	float vv_max = 0;
	// assign velocities to cell placeholders
	for (unsigned int i = 0; i < nplace; ++i) {
	    // particle number
	    const int n = h_cell.n[i];
	    if (IS_REAL_PARTICLE(n)) {
		h_cell.v[i] = make_float(h_sample.v[n]);
		// calculate maximum squared velocity
		vv_max = std::max(vv_max, h_sample.v[n] * h_sample.v[n]);
	    }
	}
	// set maximum velocity magnitude
	v_max = std::sqrt(vv_max);
	// set sum over maximum velocity magnitudes to zero
	v_max_sum = v_max;
	// copy particle velocities from host to GPU (after force calculation!)
	cuda::copy(h_cell.v, g_cell.v, stream_);
#else
	// copy periodically reduced particle positions from host to GPU
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.r[i] = make_float(h_sample.r[i]);
	}
	cuda::copy(h_part.r, g_part.r, stream_);
	// replicate to periodically extended particle positions
	cuda::copy(g_part.r, g_part.R, stream_);
	// calculate forces, potential energy and virial equation sum
	cuda::configure(dim_.grid, dim_.block, dim_.threads_per_block() * sizeof(g_vector), stream_);
	gpu::ljfluid::mdstep(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);

	// maximum squared velocity
	float vv_max = 0;
	// copy particle velocities from host to GPU (after force calculation!)
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.v[i] = make_float(h_sample.v[i]);
	    // calculate maximum squared velocity
	    vv_max = std::max(vv_max, h_sample.v[i] * h_sample.v[i]);
	}
	// set maximum velocity magnitude
	v_max = std::sqrt(vv_max);
	cuda::copy(h_part.v, g_part.v, stream_);
#endif
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to restore system state from phase space sample");
    }

#ifdef USE_CELL
    // CUDA time for cell lists initialisation
    m_times[6] += event_[1] - event_[0];
#endif
}

/**
 * seed random number generator
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::rng(unsigned int seed)
{
    LOG("random number generator seed: " << seed);
    try {
	rng_.set(seed, stream_);
	stream_.synchronize();
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

    // particles per 2- or 3-dimensional unit cell
    const unsigned int m = 2 * (dimension - 1);
    // lower boundary for number of particles per lattice dimension
    unsigned int n = std::pow(npart / m, 1.f / dimension);
    // lower boundary for total number of lattice sites
    unsigned int N = m * std::pow(n, dimension);

    if (N < npart) {
	n += 1;
	N = m * std::pow(n, dimension);
    }
    if (N > npart) {
	LOG_WARNING("lattice not fully occupied (" << N << " sites)");
    }

    // minimum distance in 2- or 3-dimensional fcc lattice
    LOG("minimum lattice distance: " << (box_ / n) / std::sqrt(2.f));

    try {
#ifdef USE_CELL
	g_cell.r.reserve(dim_.threads());
	// compute particle lattice positions on GPU
	event_[0].record(stream_);
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::ljfluid::lattice(g_cell.r, n);
	event_[1].record(stream_);
	// TODO randomly permute particles to increase force summing accuracy
	cuda::copy(g_cell.r, g_cell2.r, stream_);
	// assign particles to cells
	event_[2].record(stream_);
	cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	gpu::ljfluid::assign_cells(g_cell2.r, g_cell.r, g_cell.n);
	event_[3].record(stream_);
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0;
	// replicate particle positions to periodically extended positions
	cuda::copy(g_cell.r, g_cell.R, stream_);
	// calculate forces, potential energy and virial equation sum
	cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	gpu::ljfluid::mdstep(g_cell.r, g_cell.v, g_cell.f, g_cell.n, g_cell.en, g_cell.virial);
#else
	// compute particle lattice positions on GPU
	event_[0].record(stream_);
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::ljfluid::lattice(g_part.r, n);
	event_[1].record(stream_);
	// TODO randomly permute particles to increase force summing accuracy
	// replicate particle positions to periodically extended positions
	cuda::copy(g_part.r, g_part.R, stream_);
	// calculate forces, potential energy and virial equation sum
	cuda::configure(dim_.grid, dim_.block, dim_.threads_per_block() * sizeof(g_vector), stream_);
	gpu::ljfluid::mdstep(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
#endif

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute particle lattice positions on GPU");
    }

    // CUDA time for lattice generation
    m_times[4] += event_[1] - event_[0];
#ifdef USE_CELL
    // CUDA time for cell lists initialisation
    m_times[6] += event_[3] - event_[2];
#endif
}

/**
 * set system temperature according to Maxwell-Boltzmann distribution
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::temperature(float temp)
{
    LOG("initialising velocities from Maxwell-Boltzmann distribution at temperature: " << temp);
    try {
#ifdef USE_CELL
	g_cell.v.reserve(dim_.threads());
	// set velocities using Maxwell-Boltzmann distribution at temperature
	event_[0].record(stream_);
	rng_.boltzmann(g_cell.v, temp, stream_);
	event_[1].record(stream_);
	// copy particle velocities from GPU to host
	cuda::copy(g_cell.v, h_cell.v, stream_);
	stream_.synchronize();
	for (unsigned int i = 0; i < npart; ++i) {
	    h_sample.v[i] = T(h_cell.v[i]);
	}
	// copy particle number tags from GPU to host
	cuda::copy(g_cell.n, h_cell.n, stream_);
#else
	// set velocities using Maxwell-Boltzmann distribution at temperature
	event_[0].record(stream_);
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::ljfluid::boltzmann(g_part.v, temp, rng_.state());
	event_[1].record(stream_);
	// copy particle velocities from GPU to host
	cuda::copy(g_part.v, h_part.v, stream_);
	stream_.synchronize();
	for (unsigned int i = 0; i < npart; ++i) {
	    h_sample.v[i] = T(h_part.v[i]);
	}
#endif
	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute Maxwell-Boltzmann distributed velocities on GPU");
    }

    // CUDA time for Maxwell-Boltzmann distribution
    m_times[5] += event_[1] - event_[0];

    // compute center of mass velocity
    T v_cm = mean(h_sample.v.begin(), h_sample.v.end());
    // set center of mass velocity to zero
    for (unsigned int i = 0; i < npart; ++i) {
	h_sample.v[i] -= v_cm;
    }

    try {
#ifdef USE_CELL
	// maximum squared velocity
	float vv_max = 0;
	// assign velocities to cell placeholders
	for (unsigned int i = 0; i < nplace; ++i) {
	    // particle number
	    const int n = h_cell.n[i];
	    if (IS_REAL_PARTICLE(n)) {
		// assign velocity to cell placeholder
		h_cell.v[i] = make_float(h_sample.v[n]);
		// calculate maximum squared velocity
		vv_max = std::max(vv_max, h_sample.v[n] * h_sample.v[n]);
	    }
	}
	// set maximum velocity magnitude
	v_max = std::sqrt(vv_max);
	// initialise sum over maximum velocity magnitudes since last cell lists update
	v_max_sum = v_max;
	// copy particle velocities from host to GPU
	cuda::copy(h_cell.v, g_cell.v, stream_);
#else
	// maximum squared velocity
	float vv_max = 0;
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.v[i] = make_float(h_sample.v[i]);
	    // calculate maximum squared velocity
	    vv_max = std::max(vv_max, h_sample.v[i] * h_sample.v[i]);
	}
	// set maximum velocity magnitude
	v_max = std::sqrt(vv_max);
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
 * write parameters to HDF5 parameter group
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::attrs(H5::Group const& param) const
{
    H5xx::group node(param.createGroup("mdsim"));
    node["dimension"] = dimension;
    node["particles"] = npart;
    node["blocks"] = dim_.blocks_per_grid();
    node["threads"] = dim_.threads_per_block();
    node["density"] = density_;
    node["box_length"] = box_;
    node["timestep"] = timestep_;
    node["cutoff_distance"] = r_cut;
#ifdef USE_CELL
    node["cells"] = ncell;
    node["placeholders"] = nplace;
    node["cell_length"] = cell_length_;
    node["cell_occupancy"] = cell_occupancy_;
#endif
#ifdef USE_SMOOTH_POTENTIAL
    node["smooth_distance"] = r_smooth;
#endif
}

/**
 * stream MD simulation step on GPU
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::mdstep()
{
    event_[1].record(stream_);
#ifdef USE_CELL
    // first leapfrog step of integration of differential equations of motion
    try {
	cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	gpu::ljfluid::inteq(g_cell.r, g_cell.R, g_cell.v, g_cell.f);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream first leapfrog step on GPU");
    }
    event_[2].record(stream_);

    // update cell lists
    if (v_max_sum * timestep_ > r_skin / 2) {
	try {
	    cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	    gpu::ljfluid::update_cells(g_cell.r, g_cell.R, g_cell.v, g_cell.n, g_cell2.r, g_cell2.R, g_cell2.v, g_cell2.n);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to stream cell list update on GPU");
	}
	event_[3].record(stream_);

	try {
	    cuda::copy(g_cell2.r, g_cell.r, stream_);
	    cuda::copy(g_cell2.R, g_cell.R, stream_);
	    cuda::copy(g_cell2.v, g_cell.v, stream_);
	    cuda::copy(g_cell2.n, g_cell.n, stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to replicate cell lists on GPU");
	}
    }
    event_[4].record(stream_);

    // Lennard-Jones force calculation
    try {
	cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	gpu::ljfluid::mdstep(g_cell.r, g_cell.v, g_cell.f, g_cell.n, g_cell.en, g_cell.virial);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream force calculation on GPU");
    }
#else
    // first leapfrog step of integration of differential equations of motion
    try {
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::ljfluid::inteq(g_part.r, g_part.R, g_part.v, g_part.f);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream first leapfrog step on GPU");
    }
    event_[2].record(stream_);

    // Lennard-Jones force calculation
    try {
	cuda::configure(dim_.grid, dim_.block, dim_.threads_per_block() * sizeof(g_vector), stream_);
	gpu::ljfluid::mdstep(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream force calculation on GPU");
    }
#endif
    event_[0].record(stream_);
}

/**
 * synchronize MD simulation step on GPU
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::synchronize()
{
    try {
	// wait for MD simulation step on GPU to finish
	event_[0].synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("MD simulation step on GPU failed");
    }

    // CUDA time for MD simulation step
    m_times[0] += event_[0] - event_[1];
    // CUDA time for velocity-Verlet integration
    m_times[1] += event_[2] - event_[1];
#ifdef USE_CELL
    // CUDA time for Lennard-Jones force update
    m_times[2] += event_[0] - event_[4];

    if (v_max_sum * timestep_ > r_skin / 2) {
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0;
	// CUDA time for cell lists update
	m_times[7] += event_[3] - event_[2];
	// CUDA time for cell lists memcpy
	m_times[8] += event_[4] - event_[3];
    }
#else
    // CUDA time for Lennard-Jones force update
    m_times[2] += event_[0] - event_[2];
#endif
}

/**
 * copy MD simulation step results from GPU to host
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::sample()
{
    // mean potential energy per particle
    h_sample.en_pot = 0;
    // mean virial equation sum per particle
    h_sample.virial = 0;

#ifdef USE_CELL
    // copy MD simulation step results from GPU to host
    try {
	event_[1].record(stream_);
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
	event_[0].record(stream_);

	// wait for CUDA operations to finish
	event_[0].synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy MD simulation step results from GPU to host");
    }

    // number of particles found in cells
    unsigned int count = 0;
    // maximum squared velocity
    float vv_max = 0;

    for (unsigned int i = 0; i < nplace; ++i) {
	// particle number
	const int n = h_cell.n[i];
	// check if real particle
	if (IS_REAL_PARTICLE(n)) {
	    // copy periodically reduced particle positions
	    h_sample.r[n] = T(h_cell.r[i]);
	    // copy periodically extended particle positions
	    h_sample.R[n] = T(h_cell.R[i]);
	    // copy particle velocities
	    h_sample.v[n] = T(h_cell.v[i]);
	    // calculate mean potential energy per particle
	    h_sample.en_pot += (h_cell.en[i] - h_sample.en_pot) / ++count;
	    // calculate mean virial equation sum per particle
	    h_sample.virial += (h_cell.virial[i] - h_sample.virial) / count;
	    // calculate maximum squared velocity
	    vv_max = std::max(vv_max, h_sample.v[n] * h_sample.v[n]);
	}
    }
    // set maximum velocity magnitude
    v_max = std::sqrt(vv_max);
    // add to sum over maximum velocity magnitudes since last cell lists update
    v_max_sum += v_max;
    // validate number of particles
    if (count != npart) {
	throw exception("particle loss while updating cell lists");
    }
#else
    // copy MD simulation step results from GPU to host
    try {
	event_[1].record(stream_);
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
	event_[0].record(stream_);

	// wait for CUDA operations to finish
	event_[0].synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy MD simulation step results from GPU to host");
    }

    // maximum squared velocity
    float vv_max = 0;

    for (unsigned int n = 0; n < npart; ++n) {
	// copy periodically reduced particle positions
	h_sample.r[n] = T(h_part.r[n]);
	// copy periodically extended particle positions
	h_sample.R[n] = T(h_part.R[n]);
	// copy particle velocities
	h_sample.v[n] = T(h_part.v[n]);
	// calculate mean potential energy per particle
	h_sample.en_pot += (h_part.en[n] - h_sample.en_pot) / (n + 1);
	// calculate mean virial equation sum per particle
	h_sample.virial += (h_part.virial[n] - h_sample.virial) / (n + 1);
	// calculate maximum squared velocity
	vv_max = std::max(vv_max, h_sample.v[n] * h_sample.v[n]);
    }
    // set maximum velocity magnitude
    v_max = std::sqrt(vv_max);
#endif

    // ensure that system is still in valid state after MD step
    if (std::isnan(h_sample.en_pot)) {
	throw exception("potential energy diverged due to excessive timestep or density");
    }

    // CUDA time for sample memcpy
    m_times[3] += event_[0] - event_[1];
}

/**
 * returns and resets CUDA time statistics
 */
template <unsigned dimension, typename T>
perf_counters ljfluid<dimension, T>::times()
{
    perf_counters times(m_times);
    // reset performance counters
    for (unsigned int i = 0; i < m_times.size(); ++i) {
	m_times[i].clear();
    }
    return times;
}

} // namespace mdsim

#undef foreach

#endif /* ! MDSIM_LJFLUID_HPP */
