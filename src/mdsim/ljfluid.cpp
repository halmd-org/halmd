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

#include <algorithm>
#include <boost/foreach.hpp>
#include <cmath>
#include "exception.hpp"
#include "gpu/ljfluid_glue.hpp"
#include "ljfluid.hpp"
#include "log.hpp"
#include "statistics.hpp"

#define foreach BOOST_FOREACH

namespace mdsim
{

/**
 * initialise fixed simulation parameters
 */
ljfluid::ljfluid()
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
void ljfluid::particles(unsigned int value)
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

    // allocate global device memory for system state
    try {
	g_part.r.resize(npart);
	g_part.R.resize(npart);
	g_part.v.resize(npart);
	g_part.f.resize(npart);
	g_part.tag.resize(npart);
	g_part.en.resize(npart);
	g_part.virial.resize(npart);

	g_part.v_max.resize(REDUCE_BLOCKS);
	g_part.en_sum.resize(REDUCE_BLOCKS);
	g_part.virial_sum.resize(REDUCE_BLOCKS);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for system state");
    }

#ifdef USE_CELL
    // allocate global device memory for sorting buffers
    try {
	g_sort.r.resize(npart);
	g_sort.R.resize(npart);
	g_sort.v.resize(npart);
	g_sort.tag.resize(npart);
	g_aux.cell.resize(npart);
	g_aux.idx.resize(npart);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for sorting buffers");
    }
#endif

    // allocate page-locked host memory for system state
    try {
	h_part.r.resize(npart);
	h_part.R.resize(npart);
	h_part.v.resize(npart);
	h_part.tag.resize(npart);
	// particle forces reside only in GPU memory

	h_part.v_max.resize(REDUCE_BLOCKS);
	h_part.en_sum.resize(REDUCE_BLOCKS);
	h_part.virial_sum.resize(REDUCE_BLOCKS);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate page-locked host memory for system state");
    }
}

/**
 * set particle density
 */
void ljfluid::density(float value)
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
void ljfluid::box(float value)
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
void ljfluid::cell_occupancy(float value)
{
    LOG("desired average cell occupancy: " << value);

    // fixed cell size due to fixed number of CUDA execution threads per block
    cell_size_ = CELL_SIZE;
    LOG("number of placeholders per cell: " << cell_size_);

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
    LOG("total number of cell placeholders: " << nplace);

    // set effective average cell occupancy
    cell_occupancy_ = npart * 1. / nplace;
    LOG("effective average cell occupancy: " << cell_occupancy_);

    if (cell_occupancy_ > 1.) {
	throw exception("average cell occupancy must not be larger than 1.0");
    }
    else if (cell_occupancy_ > 0.5) {
	LOG_WARNING("average cell occupancy is larger than 0.5");
    }

    // volume of n-dimensional sphere with neighbour list radius
    const float nbl_sphere = ((dimension + 1) * M_PI / 3) * std::pow(cell_length_, dimension);
    // set number of placeholders per neighbour list
    nbl_size = std::ceil((density_ / value) * nbl_sphere);
    LOG("number of placeholders per neighbour list: " << nbl_size);

    // copy cell parameters to device symbols
    try {
	cuda::copy(ncell, gpu::ljfluid::ncell);
	cuda::copy(cell_length_, gpu::ljfluid::r_cell);
	cuda::copy(std::pow(cell_length_, 2), gpu::ljfluid::rr_cell);
	cuda::copy(nbl_size, gpu::ljfluid::nbl_size);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy cell parameters to device symbols");
    }

    // set Hilbert space-filling curve recursion level
#ifdef DIM_3D
    const uint level = std::min(10.f, ceilf(logf(box_) / M_LN2));
#else
    const uint level = std::min(16.f, ceilf(logf(box_) / M_LN2));
#endif
    LOG("Hilbert space-filling curve recursion level: " << level);
    try {
	cuda::copy(level, gpu::ljfluid::sfc_level);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy Hilbert curve recursion level to device symbol");
    }
}
#endif /* USE_CELL */

/**
 * set number of CUDA execution threads
 */
void ljfluid::threads(unsigned int value)
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

    // allocate global device memory for cell placeholders
    try {
	g_cell.resize(dim_cell_.threads());
	g_nbl.resize(npart * nbl_size);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory cell placeholders");
    }

#endif /* USE_CELL */

    // allocate global device memory for placeholder particles
    try {
	g_part.r.reserve(dim_.threads());
	g_part.R.reserve(dim_.threads());
	g_part.v.reserve(dim_.threads());
	g_part.f.reserve(dim_.threads());
	g_part.tag.reserve(dim_.threads());
	g_part.en.reserve(dim_.threads());
	g_part.virial.reserve(dim_.threads());
#ifdef USE_CELL
	g_nbl.reserve(dim_.threads() * nbl_size);
	cuda::copy(uint(dim_.threads()), gpu::ljfluid::nbl_stride);
#endif
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for placeholder particles");
    }

#ifdef USE_CELL
    // bind GPU textures to global device memory arrays
    try {
	gpu::ljfluid::r.bind(g_part.r);
	gpu::ljfluid::v.bind(g_part.v);
	gpu::ljfluid::R.bind(g_part.R);
	gpu::ljfluid::tag.bind(g_part.tag);
    }
    catch (cuda::error const& e) {
	throw exception("failed to bind GPU textures to global device memory arrays");
    }

    // allocate global device memory for sorting buffers
    try {
	g_sort.r.reserve(dim_.threads());
	g_sort.R.reserve(dim_.threads());
	g_sort.v.reserve(dim_.threads());
	g_sort.tag.reserve(dim_.threads());
	g_aux.cell.reserve(dim_.threads());
	g_aux.offset.resize(dim_cell_.blocks_per_grid());
	g_aux.idx.reserve(dim_.threads());
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for sorting buffers");
    }

    // allocate global device memory for radix sort
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
#endif

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
void ljfluid::restore(trajectory_sample::visitor visitor)
{
    // read phase space sample
    visitor(h_sample.r, h_sample.v);

    try {
	// assign particle tags
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::ljfluid::init_tags(g_part.tag);
	// copy periodically reduced particle positions from host to GPU
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.r[i] = make_float(h_sample.r[i]);
	}
	cuda::copy(h_part.r, g_part.r, stream_);
#ifdef USE_CELL
#ifdef USE_HILBERT_ORDER
	// order particles after Hilbert space-filling curve
	hilbert_order(stream_);
#endif
	// assign particles to cells
	assign_cells(stream_);
	// update neighbour lists
	update_neighbours(stream_);
#endif
	// calculate forces
	update_forces(stream_);
	// calculate potential energy
	potential_energy(stream_);
	// calculate virial equation sum
	virial_sum(stream_);
	// replicate to periodically extended particle positions
	cuda::copy(g_part.r, g_part.R, stream_);

	// copy particle velocities from host to GPU (after force calculation!)
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.v[i] = make_float(h_sample.v[i]);
	}
	cuda::copy(h_part.v, g_part.v, stream_);
#ifdef USE_CELL
	// calculate maximum velocity magnitude
	maximum_velocity(stream_);
#endif

	// wait for GPU operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to restore system state from phase space sample");
    }

    // set initial sum over maximum velocity magnitudes since last cell lists update
    v_max_sum = *std::max_element(h_part.v_max.begin(), h_part.v_max.end());
}

/**
 * seed random number generator
 */
void ljfluid::rng(unsigned int seed)
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
void ljfluid::rng(mdsim::rand48::state_type const& state)
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
void ljfluid::lattice()
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
	// assign particle tags
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::ljfluid::init_tags(g_part.tag);
	// compute particle lattice positions on GPU
	event_[0].record(stream_);
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::ljfluid::lattice(g_part.r, n);
	event_[1].record(stream_);
#ifdef USE_CELL
#ifdef USE_HILBERT_ORDER
	// order particles after Hilbert space-filling curve
	hilbert_order(stream_);
#endif
	// assign particles to cells
	assign_cells(stream_);
	// update neighbour lists
	update_neighbours(stream_);
#endif
	// calculate forces
	update_forces(stream_);
	// calculate potential energy
	potential_energy(stream_);
	// calculate virial equation sum
	virial_sum(stream_);
	// replicate particle positions to periodically extended positions
	cuda::copy(g_part.r, g_part.R, stream_);

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute particle lattice positions on GPU");
    }

    // GPU time for lattice generation
    m_times[GPU_TIME_LATTICE] += event_[1] - event_[0];
}

/**
 * set system temperature according to Maxwell-Boltzmann distribution
 */
void ljfluid::temperature(float temp)
{
    LOG("initialising velocities from Maxwell-Boltzmann distribution at temperature: " << temp);
    try {
	// set velocities using Maxwell-Boltzmann distribution at temperature
	event_[0].record(stream_);
	rng_.boltzmann(g_part.v, temp, stream_);
	event_[1].record(stream_);
	// copy particle velocities from GPU to host
	cuda::copy(g_part.v, h_part.v, stream_);
	stream_.synchronize();
	for (unsigned int i = 0; i < npart; ++i) {
	    h_sample.v[i] = h_part.v[i];
	}
	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute Maxwell-Boltzmann distributed velocities on GPU");
    }

    // GPU time for Maxwell-Boltzmann distribution
    m_times[GPU_TIME_BOLTZMANN] += event_[1] - event_[0];

    // compute center of mass velocity
    hvector v_cm = mean(h_sample.v.begin(), h_sample.v.end());
    // set center of mass velocity to zero
    for (unsigned int i = 0; i < npart; ++i) {
	h_sample.v[i] -= v_cm;
    }

    try {
	// copy particle velocities from host to GPU
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.v[i] = make_float(h_sample.v[i]);
	}
	cuda::copy(h_part.v, g_part.v, stream_);
#ifdef USE_CELL
	// calculate maximum velocity magnitude
	maximum_velocity(stream_);
#endif

	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to set center of mass velocity to zero");
    }

    // set initial sum over maximum velocity magnitudes since last cell lists update
    v_max_sum = *std::max_element(h_part.v_max.begin(), h_part.v_max.end());
}

/**
 * set simulation timestep
 */
void ljfluid::timestep(float value)
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
void ljfluid::attrs(H5::Group const& param) const
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
    node["neighbours"] = nbl_size;
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
void ljfluid::mdstep()
{
    event_[1].record(stream_);

    // first leapfrog step of integration of differential equations of motion
    try {
	velocity_verlet(stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream first leapfrog step on GPU");
    }
    event_[2].record(stream_);

#ifdef USE_CELL
    // update cell lists
    if (v_max_sum * timestep_ > r_skin / 2) {
#ifdef USE_HILBERT_ORDER
	try {
	    hilbert_order(stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to stream hilbert space-filling curve sort on GPU");
	}
#endif
	event_[5].record(stream_);

	try {
	    assign_cells(stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to stream cell list update on GPU");
	}
	event_[6].record(stream_);

	try {
	    update_neighbours(stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to stream neighbour lists update on GPU");
	}
    }
    event_[7].record(stream_);
#endif

    // Lennard-Jones force calculation
    try {
	update_forces(stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream force calculation on GPU");
    }
    event_[3].record(stream_);

    // potential energy sum calculation
    try {
	potential_energy(stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream potential energy sum calculation on GPU");
    }
    event_[4].record(stream_);

    // virial equation sum calculation
    try {
	virial_sum(stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream virial equation sum calculation on GPU");
    }

#ifdef USE_CELL
    event_[8].record(stream_);

    // maximum velocity calculation
    try {
	maximum_velocity(stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream maximum velocity calculation on GPU");
    }
#endif
    event_[0].record(stream_);
}

/**
 * synchronize MD simulation step on GPU
 */
void ljfluid::synchronize()
{
    try {
	// wait for MD simulation step on GPU to finish
	event_[0].synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("MD simulation step on GPU failed");
    }

    // GPU time for MD simulation step
    m_times[GPU_TIME_MDSTEP] += event_[0] - event_[1];
    // GPU time for velocity-Verlet integration
    m_times[GPU_TIME_VELOCITY_VERLET] += event_[2] - event_[1];
#ifdef USE_CELL
    // GPU time for Lennard-Jones force update
    m_times[GPU_TIME_UPDATE_FORCES] += event_[3] - event_[7];
    // GPU time for potential energy sum calculation
    m_times[GPU_TIME_POTENTIAL_ENERGY] += event_[4] - event_[3];
    // GPU time for virial equation sum calculation
    m_times[GPU_TIME_VIRIAL_SUM] += event_[8] - event_[4];
    // GPU time for maximum velocity calculation
    m_times[GPU_TIME_MAXIMUM_VELOCITY] += event_[0] - event_[8];

    if (v_max_sum * timestep_ > r_skin / 2) {
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0;
#ifdef USE_HILBERT_ORDER
	// GPU time for Hilbert curve sort
	m_times[GPU_TIME_HILBERT_SORT] += event_[5] - event_[2];
#endif
	// GPU time for cell lists update
	m_times[GPU_TIME_UPDATE_CELLS] += event_[6] - event_[5];
	// GPU time for neighbour lists update
	m_times[GPU_TIME_UPDATE_NEIGHBOURS] += event_[7] - event_[6];
    }
#else
    // GPU time for Lennard-Jones force update
    m_times[GPU_TIME_UPDATE_FORCES] += event_[3] - event_[2];
    // GPU time for potential energy sum calculation
    m_times[GPU_TIME_POTENTIAL_ENERGY] += event_[4] - event_[3];
    // GPU time for virial equation sum calculation
    m_times[GPU_TIME_VIRIAL_SUM] += event_[0] - event_[4];
#endif

    // mean potential energy per particle
    h_sample.en_pot = 0;
    for (unsigned int i = 0; i < h_part.en_sum.size(); ++i) {
	h_sample.en_pot += (double) h_part.en_sum[i].x + (double) h_part.en_sum[i].y;
    }
    h_sample.en_pot /= npart;
    // mean virial equation sum per particle
    h_sample.virial = 0;
    for (unsigned int i = 0; i < h_part.virial_sum.size(); ++i) {
	h_sample.virial += (double) h_part.virial_sum[i].x + (double) h_part.virial_sum[i].y;
    }
    h_sample.virial /= npart;

    // ensure that system is still in valid state after MD step
    if (!std::isfinite(h_sample.en_pot)) {
	throw exception("potential energy diverged");
    }

#ifdef USE_CELL
    // add to sum over maximum velocity magnitudes since last cell lists update
    v_max_sum += *std::max_element(h_part.v_max.begin(), h_part.v_max.end());
#endif
}

/**
 * copy MD simulation step results from GPU to host
 */
void ljfluid::sample()
{
    // copy MD simulation step results from GPU to host
    try {
	event_[1].record(stream_);
	// copy periodically reduce particles positions
	cuda::copy(g_part.r, h_part.r, stream_);
	// copy periodically extended particles positions
	cuda::copy(g_part.R, h_part.R, stream_);
	// copy particle velocities
	cuda::copy(g_part.v, h_part.v, stream_);
	// copy particle tags
	cuda::copy(g_part.tag, h_part.tag, stream_);
	event_[0].record(stream_);

	// wait for CUDA operations to finish
	event_[0].synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy MD simulation step results from GPU to host");
    }

    for (unsigned int j = 0; j < npart; ++j) {
	// particle tracking number
	const int n = h_part.tag[j];
	// copy periodically reduced particle positions
	h_sample.r[n] = h_part.r[j];
	// copy periodically extended particle positions
	h_sample.R[n] = h_part.R[j];
	// copy particle velocities
	h_sample.v[n] = h_part.v[j];
    }

    // GPU time for sample memcpy
    m_times[GPU_TIME_SAMPLE_MEMCPY] += event_[0] - event_[1];
}

/**
 * returns and resets GPU time statistics
 */
perf_counters ljfluid::times()
{
    perf_counters times(m_times);
    // reset performance counters
    for (unsigned int i = 0; i < m_times.size(); ++i) {
	m_times[i].clear();
    }
    return times;
}

/*
 * first leapfrog step of integration of differential equations of motion
 */
void ljfluid::velocity_verlet(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid::inteq(g_part.r, g_part.R, g_part.v, g_part.f);
}

/*
 * Lennard-Jones force calculation
 */
void ljfluid::update_forces(cuda::stream& stream)
{
#ifdef USE_CELL
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid::mdstep(g_part.r, g_part.v, g_part.f, g_nbl, g_part.en, g_part.virial);
#else
    cuda::configure(dim_.grid, dim_.block, dim_.threads_per_block() * sizeof(gvector), stream);
    gpu::ljfluid::mdstep(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
#endif
}

/*
 * potential energy sum calculation
 */
void ljfluid::potential_energy(cuda::stream& stream)
{
    cuda::configure(REDUCE_BLOCKS, REDUCE_THREADS, REDUCE_THREADS * sizeof(float2), stream_);
    gpu::ljfluid::potential_energy_sum(g_part.en, g_part.en_sum);
    cuda::copy(g_part.en_sum, h_part.en_sum, stream_);
}

/*
 * virial equation sum calculation
 */
void ljfluid::virial_sum(cuda::stream& stream)
{
    cuda::configure(REDUCE_BLOCKS, REDUCE_THREADS, REDUCE_THREADS * sizeof(float2), stream);
    gpu::ljfluid::potential_energy_sum(g_part.virial, g_part.virial_sum);
    cuda::copy(g_part.virial_sum, h_part.virial_sum, stream);
}

#ifdef USE_CELL

/*
 * maximum velocity calculation
 */
void ljfluid::maximum_velocity(cuda::stream& stream)
{
    cuda::configure(REDUCE_BLOCKS, REDUCE_THREADS, REDUCE_THREADS * sizeof(float), stream);
    gpu::ljfluid::maximum_velocity(g_part.v, g_part.v_max);
    cuda::copy(g_part.v_max, h_part.v_max, stream);
}

/**
 * assign particles to cells
 */
void ljfluid::assign_cells(cuda::stream& stream)
{
    // compute cell indices for particle positions
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid::compute_cell(g_part.r, g_aux.cell);

    // generate permutation array
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid::gen_index(g_aux.idx);
    radix_sort(g_aux.cell, g_aux.idx, stream);

    // compute global cell offsets in sorted particle list
    cuda::memset(g_aux.offset, 0xff);
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid::find_cell_offset(g_aux.cell, g_aux.offset);

    // assign particles to cells
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream);
    gpu::ljfluid::assign_cells(g_aux.cell, g_aux.offset, g_aux.idx, g_cell);
}

/*
 * update neighbour lists
 */
void ljfluid::update_neighbours(cuda::stream& stream)
{
    // mark neighbour list placeholders as virtual particles
    cuda::memset(g_nbl, 0xff);
    // build neighbour lists
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
    gpu::ljfluid::update_neighbours(g_cell, g_nbl);
}

#ifdef USE_HILBERT_ORDER

/**
 * order particles after Hilbert space-filling curve
 */
void ljfluid::hilbert_order(cuda::stream& stream)
{
    // compute Hilbert space-filling curve for particles
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid::sfc_hilbert_encode(g_part.r, g_aux.cell);

    // generate permutation array
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid::gen_index(g_aux.idx);
    radix_sort(g_aux.cell, g_aux.idx, stream);

    // order particles by permutation
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid::order_particles(g_aux.idx, g_sort.r, g_sort.R, g_sort.v, g_sort.tag);
    cuda::copy(g_sort.r, g_part.r, stream);
    cuda::copy(g_sort.R, g_part.R, stream);
    cuda::copy(g_sort.v, g_part.v, stream);
    cuda::copy(g_sort.tag, g_part.tag, stream);
}

#endif /* USE_HILBERT_ORDER */

#endif /* USE_CELL */

} // namespace mdsim
