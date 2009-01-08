/* Lennard-Jones fluid simulation using CUDA
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

#ifndef LJGPU_LJFLUID_LJFLUID_GPU_NBR_HPP
#define LJGPU_LJFLUID_LJFLUID_GPU_NBR_HPP

#include <algorithm>
#include <boost/foreach.hpp>
#include <cuda_wrapper.hpp>
#include <ljgpu/algorithm/radix_sort.hpp>
#include <ljgpu/algorithm/reduce.hpp>
#include <ljgpu/math/stat.hpp>
#include <ljgpu/ljfluid/gpu/hilbert.hpp>
#include <ljgpu/ljfluid/gpu/lattice.hpp>
#include <ljgpu/ljfluid/gpu/ljfluid_nbr.hpp>
#include <ljgpu/ljfluid/ljfluid_traits.hpp>
#include <ljgpu/rng/rand48.hpp>
#include <ljgpu/sample/perf.hpp>
#include <ljgpu/sample/sample.hpp>
#include <ljgpu/util/H5param.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <ljgpu/util/exception.hpp>
#include <ljgpu/util/log.hpp>

namespace ljgpu
{

template <int dimension>
class ljfluid_gpu_impl_neighbour
{
public:
    typedef typename ljfluid_gpu_traits<dimension>::float_type float_type;
    typedef typename ljfluid_gpu_traits<dimension>::vector_type vector_type;
    typedef typename ljfluid_gpu_traits<dimension>::gpu_vector_type gpu_vector_type;
    typedef typename ljfluid_gpu_traits<dimension>::trajectory_sample trajectory_sample;
    typedef typename trajectory_sample::visitor trajectory_visitor;

public:
    /** set desired average cell occupancy */
    void cell_occupancy(float_type value);
    /** restore system state from phase space sample */
    void restore(trajectory_visitor visitor);
    /** place particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(float_type temp);

    /** stream MD simulation step on GPU */
    void mdstep();
    /** synchronize MD simulation step on GPU */
    void synchronize();
    /** copy MD simulation step results from GPU to host */
    void sample();

    /** get number of cells per dimension */
    unsigned int cells() const { return ncell; }
    /** get total number of cell placeholders */
    unsigned int placeholders() const { return nplace; }
    /** get cell length */
    float_type cell_length() const { return cell_length_; }
    /** get effective average cell occupancy */
    float_type cell_occupancy() const { return cell_occupancy_; }
    /** get number of placeholders per cell */
    unsigned int cell_size() const { return cell_size_; }
    /** get total number of placeholders per neighbour list */
    unsigned int neighbours() const { return nbl_size; }
    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;

protected:
    /** set number of particles in system */
    void particles(unsigned int value);
    /** set potential cutoff radius */
    void cutoff_radius(float_type value);
#ifdef USE_POTENTIAL_SMOOTHING
    /** set potential smoothing function scale parameter */
    void potential_smoothing(float_type value);
#endif
    /** set number of CUDA execution threads */
    void threads();
    /* set periodic box length */
    void box(float_type value);
    /** set simulation timestep */
    void timestep(float_type value);

private:
    /** first leapfrog step of integration of differential equations of motion */
    void velocity_verlet(cuda::stream& stream);
    /** Lennard-Jones force calculation */
    void update_forces(cuda::stream& stream);
    /** assign particles to cells */
    void assign_cells(cuda::stream& stream);
    /** update neighbour lists */
    void update_neighbours(cuda::stream& stream);
#if defined(USE_HILBERT_ORDER)
    /** order particles after Hilbert space-filling curve */
    void hilbert_order(cuda::stream& stream);
#endif

protected:
    /** number of particles in system */
    unsigned int npart;
    /** particle density */
    float_type density_;
    /** periodic box length */
    float_type box_;
    /** simulation timestep */
    float_type timestep_;

    /** cutoff radius for shifted Lennard-Jones potential */
    float_type r_cut;
#ifdef USE_POTENTIAL_SMOOTHING
    /** potential smoothing function scale parameter */
    float_type r_smooth;
    /** squared inverse potential smoothing function scale parameter */
    float_type rri_smooth;
#endif
    /** squared cutoff radius */
    float_type rr_cut;
    /** potential energy at cutoff radius */
    float_type en_cut;

    /** GPU random number generator */
    rand48 rng_;

    /** CUDA execution dimensions */
    cuda::config dim_;
    /** CUDA asynchronous execution */
    cuda::stream stream_;
    /** CUDA events for kernel timing */
    boost::array<cuda::event, 9> event_;

    /** trajectory sample in swappable host memory */
    trajectory_sample m_sample;
    /** GPU time accumulators */
    perf_counters m_times;

private:
    /** CUDA execution dimensions for cell-specific kernels */
    cuda::config dim_cell_;

    /** GPU radix sort */
    radix_sort<int> radix_;
    /** potential energy sum */
    reduce<tag::sum, dfloat> reduce_en;
    /** virial equation sum */
    reduce<tag::sum, dfloat> reduce_virial;
    /** maximum absolute velocity */
    reduce<tag::max, float> reduce_v_max;

    /** number of cells per dimension */
    unsigned int ncell;
    /** total number of cell placeholders */
    unsigned int nplace;
    /** cell length */
    float_type cell_length_;
    /** effective average cell occupancy */
    float_type cell_occupancy_;
    /** number of placeholders per cell */
    unsigned int cell_size_;
    /** number of placeholders per neighbour list */
    unsigned int nbl_size;

    /** cell skin */
    float_type r_skin;
    /** sum over maximum velocity magnitudes since last cell lists update */
    float_type v_max_sum;

    /** system state in page-locked host memory */
    struct {
	/** periodically reduced particle positions */
	cuda::host::vector<gpu_vector_type> r;
	/** periodically extended particle positions */
	cuda::host::vector<gpu_vector_type> R;
	/** particle velocities */
	cuda::host::vector<gpu_vector_type> v;
	/** particle tags */
	cuda::host::vector<int> tag;
    } h_part;

    /** system state in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<gpu_vector_type> r;
	/** periodically extended particle positions */
	cuda::vector<gpu_vector_type> R;
	/** particle velocities */
	cuda::vector<gpu_vector_type> v;
	/** particle forces */
	cuda::vector<gpu_vector_type> f;
	/** particle tags */
	cuda::vector<int> tag;
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
    } g_part;

    /** double buffers for particle sorting */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<gpu_vector_type> r;
	/** periodically extended particle positions */
	cuda::vector<gpu_vector_type> R;
	/** particle velocities */
	cuda::vector<gpu_vector_type> v;
	/** particle tags */
	cuda::vector<int> tag;
    } g_sort;

    /** auxiliary device memory arrays for particle sorting */
    struct {
	/** particle cells */
	cuda::vector<uint> cell;
	/** cell offsets in sorted particle list */
	cuda::vector<int> offset;
	/** permutation indices */
	cuda::vector<int> idx;
    } g_aux;

    /** cell lists in global device memory */
    cuda::vector<int> g_cell;
    /** neighbour lists in global device memory */
    cuda::vector<int> g_nbl;
};

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::particles(unsigned int value)
{
    // copy particle number to device symbol
    try {
	cuda::copy(npart, gpu::ljfluid_neighbour::npart);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy particle number to device symbol");
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
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for system state");
    }

    // allocate global device memory for sorting buffers
    try {
	g_sort.r.resize(npart);
	g_sort.R.resize(npart);
	g_sort.v.resize(npart);
	g_sort.tag.resize(npart);
	g_aux.cell.resize(npart);
	g_aux.idx.resize(npart);
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for sorting buffers");
    }

    // allocate page-locked host memory for system state
    try {
	h_part.r.resize(npart);
	h_part.R.resize(npart);
	h_part.v.resize(npart);
	h_part.tag.resize(npart);
	// particle forces reside only in GPU memory
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate page-locked host memory for system state");
    }
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::cutoff_radius(float_type value)
{
    try {
	cuda::copy(r_cut, gpu::ljfluid_neighbour::r_cut);
	cuda::copy(rr_cut, gpu::ljfluid_neighbour::rr_cut);
	cuda::copy(en_cut, gpu::ljfluid_neighbour::en_cut);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy potential cutoff symbols");
    }
}

#ifdef USE_POTENTIAL_SMOOTHING
template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::potential_smoothing(float_type value)
{
    try {
	cuda::copy(rri_smooth, gpu::ljfluid_neighbour::rri_smooth);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy potential smoothing function scale symbol");
    }
}
#endif /* USE_POTENTIAL_SMOOTHING */

/**
 * set desired average cell occupancy
 */
template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::cell_occupancy(float_type value)
{
    LOG("desired average cell occupancy: " << value);

    // fixed cell size due to fixed number of CUDA execution threads per block
    cell_size_ = gpu::ljfluid_neighbour::CELL_SIZE;
    LOG("number of placeholders per cell: " << cell_size_);

    // optimal number of cells with given cell occupancy as upper boundary
    ncell = std::ceil(std::pow(npart / (value * cell_size_), 1.f / dimension));

    // set number of cells per dimension, respecting cutoff radius
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
    float_type nbl_sphere = ((dimension + 1) * M_PI / 3) * std::pow(cell_length_, dimension);
    // set number of placeholders per neighbour list
    nbl_size = std::ceil((density_ / value) * nbl_sphere);
    LOG("number of placeholders per neighbour list: " << nbl_size);

    // copy cell parameters to device symbols
    try {
	//cuda::copy(ncell, gpu::ljfluid_gpu_cell::ncell);
	cuda::copy(ncell, gpu::ljfluid_neighbour::ncell);
	cuda::copy(cell_length_, gpu::ljfluid_neighbour::r_cell);
	cuda::copy(std::pow(cell_length_, 2), gpu::ljfluid_neighbour::rr_cell);
	cuda::copy(nbl_size, gpu::ljfluid_neighbour::nbl_size);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy cell parameters to device symbols");
    }

#if defined(USE_HILBERT_ORDER)
    // set Hilbert space-filling curve recursion depth
    uint depth = std::min((dimension == 3) ? 10.f : 16.f, ceilf(logf(box_) / M_LN2));
    LOG("Hilbert space-filling curve recursion depth: " << depth);
    try {
	cuda::copy(box_, gpu::hilbert::box);
	cuda::copy(depth, gpu::hilbert::depth);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy Hilbert curve recursion depth to device symbol");
    }
#endif
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::threads()
{
    // set CUDA execution dimensions for cell-specific kernels
    dim_cell_ = cuda::config(dim3(powf(ncell, dimension)), dim3(cell_size_));
    LOG("number of cell CUDA execution blocks: " << dim_cell_.blocks_per_grid());
    LOG("number of cell CUDA execution threads: " << dim_cell_.threads_per_block());

    // allocate global device memory for cell placeholders
    try {
	g_cell.resize(dim_cell_.threads());
	g_nbl.resize(npart * nbl_size);
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory cell placeholders");
    }

    // allocate global device memory for placeholder particles
    try {
	g_part.r.reserve(dim_.threads());
	g_part.R.reserve(dim_.threads());
	g_part.v.reserve(dim_.threads());
	g_part.f.reserve(dim_.threads());
	g_part.tag.reserve(dim_.threads());
	g_part.en.reserve(dim_.threads());
	g_part.virial.reserve(dim_.threads());
	g_nbl.reserve(dim_.threads() * nbl_size);
	cuda::copy(uint(dim_.threads()), gpu::ljfluid_neighbour::nbl_stride);
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for placeholder particles");
    }

    // bind GPU textures to global device memory arrays
    try {
	gpu::ljfluid_neighbour::texref<gpu_vector_type>::r.bind(g_part.r);
	gpu::ljfluid_neighbour::texref<gpu_vector_type>::v.bind(g_part.v);
	gpu::ljfluid_neighbour::texref<gpu_vector_type>::R.bind(g_part.R);
	gpu::ljfluid_neighbour::texref<gpu_vector_type>::tag.bind(g_part.tag);
    }
    catch (cuda::error const&) {
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
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for sorting buffers");
    }

    // allocate global device memory for radix sort
    try {
	// compute optimal number of blocks for GeForce 8800 with 16 multiprocessors
	uint threads = dim_.threads_per_block();
	uint max_blocks = (16 * 512) / (threads * gpu::radix_sort::BUCKETS_PER_THREAD / 2);
	uint blocks = std::min((npart + 2 * threads - 1) / (2 * threads), max_blocks);

	LOG("number of CUDA blocks for radix sort: " << blocks);
	LOG("number of CUDA threads for radix sort: " << threads);

	// allocate global device memory for radix sort
	radix_.resize(npart, blocks, threads);
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for radix sort");
    }
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::box(float_type value)
{
    try {
	cuda::copy(box_, gpu::ljfluid_neighbour::box);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy periodic box length to device symbol");
    }
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::timestep(float_type value)
{
    try {
	cuda::copy(timestep_, gpu::ljfluid_neighbour::timestep);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy simulation timestep to device symbol");
    }
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::restore(trajectory_visitor visitor)
{
    // read phase space sample
    visitor(m_sample.r, m_sample.v);

    try {
	// assign particle tags
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::ljfluid_neighbour::init_tags(g_part.tag);
	// copy periodically reduced particle positions from host to GPU
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.r[i] = make_float(m_sample.r[i]);
	}
	cuda::copy(h_part.r, g_part.r, stream_);
	// replicate to periodically extended particle positions
	cuda::copy(g_part.r, g_part.R, stream_);
#ifdef USE_HILBERT_ORDER
	// order particles after Hilbert space-filling curve
	hilbert_order(stream_);
#endif
	// assign particles to cells
	assign_cells(stream_);
	// update neighbour lists
	update_neighbours(stream_);
	// calculate forces
	update_forces(stream_);
	// calculate potential energy
	reduce_en(g_part.en, stream_);
	// calculate virial equation sum
	reduce_virial(g_part.virial, stream_);

	// copy particle velocities from host to GPU (after force calculation!)
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.v[i] = make_float(m_sample.v[i]);
	}
	cuda::copy(h_part.v, g_part.v, stream_);
	// calculate maximum velocity magnitude
	reduce_v_max(g_part.v, stream_);

	// wait for GPU operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to restore system state from phase space sample");
    }

    // set initial sum over maximum velocity magnitudes since last cell lists update
    v_max_sum = reduce_v_max.value();
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::lattice()
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
	gpu::ljfluid_neighbour::init_tags(g_part.tag);
	// compute particle lattice positions on GPU
	event_[0].record(stream_);
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::lattice::fcc(g_part.r, n, box_);
	event_[1].record(stream_);
#ifdef USE_HILBERT_ORDER
	// order particles after Hilbert space-filling curve
	hilbert_order(stream_);
#endif
	// assign particles to cells
	assign_cells(stream_);
	// update neighbour lists
	update_neighbours(stream_);
	// calculate forces
	update_forces(stream_);
	// calculate potential energy
	reduce_en(g_part.en, stream_);
	// calculate virial equation sum
	reduce_virial(g_part.virial, stream_);
	// replicate particle positions to periodically extended positions
	cuda::copy(g_part.r, g_part.R, stream_);

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute particle lattice positions on GPU");
    }

    // GPU time for lattice generation
    m_times["lattice"] += event_[1] - event_[0];
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::temperature(float_type temp)
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
	    m_sample.v[i] = h_part.v[i];
	}
	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute Maxwell-Boltzmann distributed velocities on GPU");
    }

    // GPU time for Maxwell-Boltzmann distribution
    m_times["boltzmann"] += event_[1] - event_[0];

    //
    // The particle velocities need to fullfill two constraints:
    //
    //  1. center of mass velocity shall be zero
    //  2. temperature of the distribution shall equal exactly the given value
    //
    // We choose the above order because shifting the center of mass velocity
    // means altering the first moment of the velocity distribution, which in
    // consequence affects the second moment, i.e. the temperature.
    //

    // compute center of mass velocity
    vector_type v_cm = mean(m_sample.v.begin(), m_sample.v.end());
    // set center of mass velocity to zero
    for (unsigned int i = 0; i < npart; ++i) {
	m_sample.v[i] -= v_cm;
    }

    // compute mean squared velocity
    double vv = 0;
    for (unsigned int i = 0; i < npart; ++i) {
	vv += (m_sample.v[i] * m_sample.v[i] - vv) / (i + 1);
    }
    // rescale velocities to accurate temperature
    double s = sqrt(temp * dimension / vv);
    for (unsigned int i = 0; i < npart; ++i) {
	m_sample.v[i] *= s;
    }

    try {
	// copy particle velocities from host to GPU
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.v[i] = make_float(m_sample.v[i]);
	}
	cuda::copy(h_part.v, g_part.v, stream_);
	// calculate maximum velocity magnitude
	reduce_v_max(g_part.v, stream_);

	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to set center of mass velocity to zero");
    }

    // set initial sum over maximum velocity magnitudes since last cell lists update
    v_max_sum = reduce_v_max.value();
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::mdstep()
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
	event_[3].record(stream_);

	try {
	    assign_cells(stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to stream cell list update on GPU");
	}
	event_[4].record(stream_);

	try {
	    update_neighbours(stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to stream neighbour lists update on GPU");
	}
    }
    event_[5].record(stream_);

    // Lennard-Jones force calculation
    try {
	update_forces(stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream force calculation on GPU");
    }
    event_[6].record(stream_);

    // potential energy sum calculation
    try {
	reduce_en(g_part.en, stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream potential energy sum calculation on GPU");
    }
    event_[7].record(stream_);

    // virial equation sum calculation
    try {
	reduce_virial(g_part.virial, stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream virial equation sum calculation on GPU");
    }
    event_[8].record(stream_);

    // maximum velocity calculation
    try {
	reduce_v_max(g_part.v, stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream maximum velocity calculation on GPU");
    }
    event_[0].record(stream_);
}

/**
 * synchronize MD simulation step on GPU
 */
template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::synchronize()
{
    try {
	// wait for MD simulation step on GPU to finish
	event_[8].synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("MD simulation step on GPU failed");
    }

    // CUDA time for MD simulation step
    m_times["mdstep"] += event_[0] - event_[1];
    // GPU time for velocity-Verlet integration
    m_times["velocity_verlet"] += event_[2] - event_[1];
    // GPU time for Lennard-Jones force update
    m_times["update_forces"] += event_[6] - event_[5];
    // GPU time for potential energy sum calculation
    m_times["potential_energy"] += event_[7] - event_[6];
    // GPU time for virial equation sum calculation
    m_times["virial_sum"] += event_[8] - event_[7];
    // GPU time for maximum velocity calculation
    m_times["maximum_velocity"] += event_[0] - event_[8];

    if (v_max_sum * timestep_ > r_skin / 2) {
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0;
#ifdef USE_HILBERT_ORDER
	// GPU time for Hilbert curve sort
	m_times["hilbert_sort"] += event_[3] - event_[2];
#endif
	// GPU time for cell lists update
	m_times["update_cells"] += event_[4] - event_[3];
	// GPU time for neighbour lists update
	m_times["update_neighbours"] += event_[5] - event_[4];
    }

    if (!std::isfinite((double) reduce_en.value())) {
	throw exception("potential energy diverged");
    }

    // add to sum over maximum velocity magnitudes since last cell lists update
    v_max_sum += reduce_v_max.value();
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::sample()
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
	m_sample.r[n] = h_part.r[j];
	// copy periodically extended particle positions
	m_sample.R[n] = h_part.R[j];
	// copy particle velocities
	m_sample.v[n] = h_part.v[j];
    }

    // mean potential energy per particle
    m_sample.en_pot = (double) reduce_en.value() / npart;
    // mean virial equation sum per particle
    m_sample.virial = (double) reduce_virial.value() / npart;

    // GPU time for sample memcpy
    m_times["sample_memcpy"] += event_[0] - event_[1];
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::attrs(H5::Group const& param) const
{
    H5xx::group node(param.openGroup("mdsim"));
    node["cells"] = ncell;
    node["placeholders"] = nplace;
    node["neighbours"] = nbl_size;
    node["cell_length"] = cell_length_;
    node["cell_occupancy"] = cell_occupancy_;
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::velocity_verlet(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid_neighbour::inteq(g_part.r, g_part.R, g_part.v, g_part.f);
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::update_forces(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid_neighbour::mdstep(g_part.r, g_part.v, g_part.f, g_nbl, g_part.en, g_part.virial);
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::assign_cells(cuda::stream& stream)
{
    // compute cell indices for particle positions
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid_neighbour::compute_cell(g_part.r, g_aux.cell);

    // generate permutation array
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid_neighbour::gen_index(g_aux.idx);
    radix_(g_aux.cell, g_aux.idx, stream);

    // compute global cell offsets in sorted particle list
    cuda::memset(g_aux.offset, 0xff);
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid_neighbour::find_cell_offset(g_aux.cell, g_aux.offset);

    // assign particles to cells
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream);
    gpu::ljfluid_neighbour::assign_cells(g_aux.cell, g_aux.offset, g_aux.idx, g_cell);
}

template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::update_neighbours(cuda::stream& stream)
{
    // mark neighbour list placeholders as virtual particles
    cuda::memset(g_nbl, 0xff);
    // build neighbour lists
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream);
    gpu::ljfluid_neighbour::update_neighbours(g_cell, g_nbl, g_part.r);
}

#if defined(USE_HILBERT_ORDER)
template <int dimension>
void ljfluid_gpu_impl_neighbour<dimension>::hilbert_order(cuda::stream& stream)
{
    // compute Hilbert space-filling curve for particles
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::hilbert::curve(g_part.r, g_aux.cell);

    // generate permutation array
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid_neighbour::gen_index(g_aux.idx);
    radix_(g_aux.cell, g_aux.idx, stream);

    // order particles by permutation
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid_neighbour::order_particles(g_aux.idx, g_sort.r, g_sort.R, g_sort.v, g_sort.tag);
    cuda::copy(g_sort.r, g_part.r, stream);
    cuda::copy(g_sort.R, g_part.R, stream);
    cuda::copy(g_sort.v, g_part.v, stream);
    cuda::copy(g_sort.tag, g_part.tag, stream);
}
#endif /* USE_HILBERT_ORDER */

} // namespace ljgpu

#endif /* ! LJGPU_LJFLUID_LJFLUID_GPU_NBR_HPP */
