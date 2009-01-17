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

#ifndef LJGPU_MDSIM_LJFLUID_GPU_NBR_HPP
#define LJGPU_MDSIM_LJFLUID_GPU_NBR_HPP

#include <algorithm>
#include <ljgpu/algorithm/radix_sort.hpp>
#include <ljgpu/algorithm/reduce.hpp>
#include <ljgpu/mdsim/ljfluid_gpu_base.hpp>
#include <ljgpu/mdsim/gpu/hilbert.hpp>
#include <ljgpu/mdsim/gpu/lattice.hpp>
#include <ljgpu/math/stat.hpp>

namespace ljgpu
{

template <typename ljfluid_impl>
class ljfluid;

template<int dimension>
class ljfluid<ljfluid_impl_gpu_neighbour<dimension> >
    : public ljfluid_gpu_base<ljfluid_impl_gpu_neighbour<dimension> >
{
public:
    typedef ljfluid_gpu_base<ljfluid_impl_gpu_neighbour<dimension> > _Base;
    typedef gpu::ljfluid<ljfluid_impl_gpu_neighbour<dimension> > _gpu;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::gpu_vector_type gpu_vector_type;
    typedef typename _Base::sample_type sample_type;
    typedef typename sample_type::sample_visitor sample_visitor;

public:
    /** set number of particles in system */
    void particles(unsigned int value);
    /** set number of CUDA execution threads */
    void threads(unsigned int value);
    /** set desired average cell occupancy */
    void cell_occupancy(float_type value);
    /** set neighbour list skin */
    void nbl_skin(float value);

    /** restore system state from phase space sample */
    void sample(sample_visitor visitor);
    /** place particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(float_type temp);

    /** stream MD simulation step on GPU */
    void stream();
    /** synchronize MD simulation step on GPU */
    void mdstep();
    /** copy MD simulation step results from GPU to host */
    void copy();

    /** returns number of particles */
    unsigned int particles() const { return npart; }
    /** returns trajectory sample */
    sample_type const& sample() const { return m_sample; }
    /** get number of CUDA execution threads */
    unsigned int threads() const { return dim_.threads_per_block(); }
    /** get effective average cell occupancy */
    float_type cell_occupancy() const { return cell_occupancy_; }
    /** get number of cells per dimension */
    unsigned int cells() const { return ncell; }
    /** get total number of cell placeholders */
    unsigned int placeholders() const { return nplace; }
    /** get cell length */
    float_type cell_length() const { return cell_length_; }
    /** get number of placeholders per cell */
    unsigned int cell_size() const { return cell_size_; }
    /** get total number of placeholders per neighbour list */
    unsigned int neighbours() const { return nbl_size; }
    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

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

    using _Base::dim_;
    using _Base::stream_;

    /** CUDA execution dimensions for cell-specific kernels */
    cuda::config dim_cell_;
    /** CUDA events for kernel timing */
    boost::array<cuda::event, 9> event_;

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

    /** neighbour list skin */
    float_type r_skin;
    /** cutoff distance with neighbour list skin */
    float_type r_cut_skin;
    /** number of placeholders per neighbour list */
    unsigned int nbl_size;

    /** sum over maximum velocity magnitudes since last cell lists update */
    float_type v_max_sum;

    /** system state in page-locked host memory */
    struct {
	/** periodically reduced particle positions */
	cuda::host::vector<gpu_vector_type> r;
	/** periodic box traversal vectors */
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
	/** periodic box traversal vectors */
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
	/** periodic box traversal vectors */
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
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::particles(unsigned int value)
{
    _Base::particles(value);

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

/**
 * set desired average cell occupancy
 */
template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::cell_occupancy(float_type value)
{
    LOG("desired average cell occupancy: " << value);

    // fixed cell size due to fixed number of CUDA execution threads per block
    cell_size_ = _gpu::CELL_SIZE;
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

    try {
	cuda::copy(ncell, _gpu::ncell);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy cell parameters to device symbols");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::nbl_skin(float value)
{
    r_skin = value;
    r_cut_skin = r_cut + r_skin;
    LOG("neighbour list skin: " << r_skin);

    if (r_cut_skin > cell_length_) {
	throw exception("neighbour list skin is larger than cell length");
    }

    // volume of n-dimensional sphere with neighbour list radius
    float nbl_sphere = ((dimension + 1) * M_PI / 3) * std::pow(r_cut_skin, dimension);
    // set number of placeholders per neighbour list
    nbl_size = std::ceil((density_ / value) * nbl_sphere);
    LOG("number of placeholders per neighbour list: " << nbl_size);

    try {
	cuda::copy(std::pow(r_cut_skin, 2), _gpu::rr_nbl);
	cuda::copy(nbl_size, _gpu::nbl_size);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy neighbour list parameters to device symbols");
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
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::threads(unsigned int value)
{
    _Base::threads(value);

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
	cuda::copy(uint(dim_.threads()), _gpu::nbl_stride);
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for placeholder particles");
    }

    // bind GPU textures to global device memory arrays
    try {
	_gpu::r.bind(g_part.r);
	_gpu::v.bind(g_part.v);
	_gpu::R.bind(g_part.R);
	_gpu::tag.bind(g_part.tag);
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
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::sample(sample_visitor visitor)
{
    _Base::sample(visitor);

    try {
	// assign particle tags
	cuda::configure(dim_.grid, dim_.block, stream_);
	_gpu::init_tags(g_part.tag);
	// copy periodically reduced particle positions from host to GPU
	std::copy(m_sample[0].r.begin(), m_sample[0].r.end(), h_part.r.begin());
	cuda::copy(h_part.r, g_part.r, stream_);
	// set periodic box traversal vectors to zero
	cuda::memset(g_part.R, 0);
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
	    h_part.v[i] = m_sample[0].v[i];
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
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::lattice()
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
	_gpu::init_tags(g_part.tag);
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
	// set periodic box traversal vectors to zero
	cuda::memset(g_part.R, 0);

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
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::temperature(float_type temp)
{
    LOG("initialising velocities from Maxwell-Boltzmann distribution at temperature: " << temp);
    boltzmann(g_part.v, h_part.v, temp);

    try {
	reduce_v_max(g_part.v, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to calculate maximum velocity magnitude on GPU");
    }
    v_max_sum = reduce_v_max.value();
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::stream()
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
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::mdstep()
{
    try {
	// wait for MD simulation step on GPU to finish
	event_[0].synchronize();
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
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::copy()
{
    // copy MD simulation step results from GPU to host
    try {
	event_[1].record(stream_);
	// copy periodically reduce particles positions
	cuda::copy(g_part.r, h_part.r, stream_);
	// copy periodic box traversal vectors
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
	// copy periodically extended particle positions
	m_sample[0].r[n] = h_part.r[j] + (vector_type) h_part.R[j] * box_;
	// copy particle velocities
	m_sample[0].v[n] = h_part.v[j];
    }

    // mean potential energy per particle
    m_sample.en_pot = (double) reduce_en.value() / npart;
    // mean virial equation sum per particle
    m_sample.virial = (double) reduce_virial.value() / npart;

    // GPU time for sample memcpy
    m_times["sample_memcpy"] += event_[0] - event_[1];
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::velocity_verlet(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::inteq(g_part.r, g_part.R, g_part.v, g_part.f);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::update_forces(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, stream);
    if (r_smooth > 0 && thermostat_nu > 0) {
	_gpu::mdstep_smooth_nvt(g_part.r, g_part.v, g_part.f, g_nbl, g_part.en, g_part.virial);
    }
    else if (r_smooth > 0) {
	_gpu::mdstep_smooth(g_part.r, g_part.v, g_part.f, g_nbl, g_part.en, g_part.virial);
    }
    else if (thermostat_nu > 0) {
	_gpu::mdstep_nvt(g_part.r, g_part.v, g_part.f, g_nbl, g_part.en, g_part.virial);
    }
    else {
	_gpu::mdstep(g_part.r, g_part.v, g_part.f, g_nbl, g_part.en, g_part.virial);
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::assign_cells(cuda::stream& stream)
{
    // compute cell indices for particle positions
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::compute_cell(g_part.r, g_aux.cell);

    // generate permutation array
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::gen_index(g_aux.idx);
    radix_(g_aux.cell, g_aux.idx, stream);

    // compute global cell offsets in sorted particle list
    cuda::memset(g_aux.offset, 0xff);
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::find_cell_offset(g_aux.cell, g_aux.offset);

    // assign particles to cells
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream);
    _gpu::assign_cells(g_aux.cell, g_aux.offset, g_aux.idx, g_cell);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::update_neighbours(cuda::stream& stream)
{
    // mark neighbour list placeholders as virtual particles
    cuda::memset(g_nbl, 0xff);
    // build neighbour lists
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream);
    _gpu::update_neighbours(g_cell, g_nbl, g_part.r);
}

#if defined(USE_HILBERT_ORDER)
template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::hilbert_order(cuda::stream& stream)
{
    // compute Hilbert space-filling curve for particles
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::hilbert::curve(g_part.r, g_aux.cell);

    // generate permutation array
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::gen_index(g_aux.idx);
    radix_(g_aux.cell, g_aux.idx, stream);

    // order particles by permutation
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::order_particles(g_aux.idx, g_sort.r, g_sort.R, g_sort.v, g_sort.tag);
    cuda::copy(g_sort.r, g_part.r, stream);
    cuda::copy(g_sort.R, g_part.R, stream);
    cuda::copy(g_sort.v, g_part.v, stream);
    cuda::copy(g_sort.tag, g_part.tag, stream);
}
#endif /* USE_HILBERT_ORDER */

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::param(H5param& param) const
{
    _Base::param(param);

    H5xx::group node(param["mdsim"]);
    node["cells"] = ncell;
    node["placeholders"] = nplace;
    node["neighbours"] = nbl_size;
    node["cell_length"] = cell_length_;
    node["cell_occupancy"] = cell_occupancy_;
    node["neighbour_skin"] = r_skin;
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_GPU_NBR_HPP */
