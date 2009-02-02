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
#include <boost/lexical_cast.hpp>
#include <ljgpu/algorithm/reduce.hpp>
#include <ljgpu/mdsim/ljfluid_gpu_base.hpp>
#include <ljgpu/mdsim/gpu/hilbert.hpp>

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
    typedef typename _Base::host_sample_type host_sample_type;
    typedef typename _Base::gpu_sample_type gpu_sample_type;
    typedef gpu_sample_type trajectory_sample_type;
    typedef boost::variant<host_sample_type, gpu_sample_type> trajectory_sample_variant;
    typedef typename _Base::energy_sample_type energy_sample_type;

    /** static implementation properties */
    typedef boost::true_type has_trajectory_gpu_sample;
    typedef boost::true_type has_energy_gpu_sample;
    typedef boost::true_type has_thermostat;

public:
    /** set number of particles in system */
    template <typename T>
    void particles(T const& value);
    /** set number of CUDA execution threads */
    void threads(unsigned int value);
    /** set desired average cell occupancy */
    void cell_occupancy(float_type value);
    /** set neighbour list skin */
    void nbl_skin(float value);

    /** restore system state from phase space sample */
    void state(host_sample_type& sample, float_type box);
    /** place particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(float_type temp);

    /** stream MD simulation step on GPU */
    void stream();
    /** synchronize MD simulation step on GPU */
    void mdstep();
    /** sample phase space on host */
    void sample(host_sample_type& sample) const;
    /** sample phase space on GPU */
    void sample(gpu_sample_type& sample) const;
    /** sample thermodynamic equilibrium properties */
    void sample(energy_sample_type& sample) const;

    /** returns number of particles */
    unsigned int particles() const { return npart; }
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
    /** assign particle positions */
    void assign_positions();
    /** assign particle velocities */
    void assign_velocities();
    /** first leapfrog step of integration of differential equations of motion */
    void velocity_verlet(cuda::stream& stream);
    /** Lennard-Jones force calculation */
    void update_forces(cuda::stream& stream);
    /** assign particles to cells */
    void assign_cells(cuda::stream& stream);
    /** update neighbour lists */
    void update_neighbours(cuda::stream& stream);
    /** order particles after Hilbert space-filling curve */
    void hilbert_order(cuda::stream& stream);
    /** generate permutation for phase space sampling */
    void permutation(cuda::stream& stream);

private:
    using _Base::box_;
    using _Base::density_;
    using _Base::dim_;
    using _Base::m_times;
    using _Base::mixture_;
    using _Base::mpart;
    using _Base::npart;
    using _Base::potential_;
    using _Base::r_cut;
    using _Base::stream_;
    using _Base::timestep_;
    using _Base::thermostat_steps;
    using _Base::thermostat_count;
    using _Base::thermostat_temp;

    /** CUDA execution dimensions for cell-specific kernels */
    cuda::config dim_cell_;
    /** CUDA execution dimensions for phase space sampling */
    std::vector<cuda::config> dim_sample;
    /** CUDA events for kernel timing */
    boost::array<cuda::event, 11> event_;

    /** GPU radix sort */
    radix_sort<unsigned int> radix_;

    using _Base::reduce_squared_velocity;
    using _Base::reduce_velocity;
    using _Base::reduce_en;
    using _Base::reduce_virial;
    using _Base::reduce_v_max;

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
    /** number of placeholders per neighbour list */
    unsigned int nbl_size;

    /** sum over maximum velocity magnitudes since last cell lists update */
    float_type v_max_sum;

    /** system state in page-locked host memory */
    struct {
	/** tagged periodically reduced particle positions */
	cuda::host::vector<float4> mutable r;
	/** periodic box traversal vectors */
	cuda::host::vector<gpu_vector_type> mutable R;
	/** particle velocities */
	cuda::host::vector<gpu_vector_type> mutable v;
	/** particle tags */
	cuda::host::vector<unsigned int> mutable tag;
    } h_part;

    /** system state in global device memory */
    struct {
	/** tagged periodically reduced particle positions */
	cuda::vector<float4> r;
	/** periodic box traversal vectors */
	cuda::vector<gpu_vector_type> R;
	/** particle velocities */
	cuda::vector<gpu_vector_type> v;
	/** particle forces */
	cuda::vector<gpu_vector_type> f;
	/** particle tags */
	cuda::vector<unsigned int> tag;
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
    } g_part;

    /** double buffers for particle sorting */
    struct {
	/** tagged periodically reduced particle positions */
	cuda::vector<float4> r;
	/** periodic box traversal vectors */
	cuda::vector<gpu_vector_type> R;
	/** particle velocities */
	cuda::vector<gpu_vector_type> v;
    } g_part_buf;

    /** auxiliary device memory arrays for particle sorting */
    struct {
	/** particle cells */
	cuda::vector<unsigned int> cell;
	/** cell offsets in sorted particle list */
	cuda::vector<unsigned int> offset;
	/** permutation indices */
	cuda::vector<unsigned int> index;
    } g_aux;

    /** cell lists in global device memory */
    cuda::vector<unsigned int> g_cell;
    /** neighbour lists in global device memory */
    cuda::vector<unsigned int> g_nbl;
};

template <int dimension>
template <typename T>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::particles(T const& value)
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
	g_part_buf.r.resize(npart);
	g_part_buf.R.resize(npart);
	g_part_buf.v.resize(npart);
	g_aux.cell.resize(npart);
	g_aux.index.resize(npart);
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
    float const r_cut_max = *std::max_element(r_cut.begin(), r_cut.end());

    // fixed cell size due to fixed number of CUDA execution threads per block
    cell_size_ = _gpu::CELL_SIZE;
    LOG("number of placeholders per cell: " << cell_size_);

    // optimal number of cells with given cell occupancy as upper boundary
    LOG("desired average cell occupancy: " << value);
    ncell = std::ceil(std::pow(npart / (value * cell_size_), 1.f / dimension));

    // set number of cells per dimension, respecting cutoff radius
    ncell = std::min(ncell, static_cast<unsigned int>(box_ / r_cut_max));
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
    cell_occupancy_ = npart * 1.f / nplace;
    LOG("effective average cell occupancy: " << cell_occupancy_);

    if (cell_occupancy_ > 1) {
	throw exception("average cell occupancy must not be larger than 1.0");
    }
    else if (cell_occupancy_ > 0.5f) {
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
    float const r_cut_max = *std::max_element(r_cut.begin(), r_cut.end());

    r_skin = std::min(value, cell_length_ - r_cut_max);
    if (r_skin < value) {
	LOG_WARNING("reducing neighbour list skin to fixed-size cell skin");
    }
    LOG("neighbour list skin: " << r_skin);

    float const r_nbl = r_cut_max + r_skin;
    // volume of n-dimensional sphere with neighbour list radius
    float const v_nbl = ((dimension + 1) * M_PI / 3) * std::pow(r_nbl, dimension);

    nbl_size = std::ceil(v_nbl * (density_ / cell_occupancy_));
    LOG("number of placeholders per neighbour list: " << nbl_size);

    try {
	cuda::copy(std::pow(r_nbl, 2), _gpu::rr_nbl);
	cuda::copy(nbl_size, _gpu::nbl_size);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy neighbour list parameters to device symbols");
    }

#if defined(USE_HILBERT_ORDER)
    // set Hilbert space-filling curve recursion depth
    unsigned int depth = std::min((dimension == 3) ? 10.f : 16.f, ceilf(logf(box_) / M_LN2));
    LOG("Hilbert space-filling curve recursion depth: " << depth);
    try {
	cuda::copy(box_, gpu::hilbert<dimension>::box);
	cuda::copy(depth, gpu::hilbert<dimension>::depth);
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
	cuda::copy(g_nbl.data(), _gpu::g_nbl);
	cuda::copy(static_cast<unsigned int>(dim_.threads()), _gpu::nbl_stride);
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for placeholder particles");
    }

    // bind GPU textures to global device memory arrays
    try {
	_gpu::r.bind(g_part.r);
	_gpu::v.bind(g_part.v);
	_gpu::R.bind(g_part.R);
    }
    catch (cuda::error const&) {
	throw exception("failed to bind GPU textures to global device memory arrays");
    }

    // allocate global device memory for sorting buffers
    try {
	g_part_buf.r.reserve(dim_.threads());
	g_part_buf.R.reserve(dim_.threads());
	g_part_buf.v.reserve(dim_.threads());
	g_aux.cell.reserve(dim_.threads());
	g_aux.offset.resize(dim_cell_.blocks_per_grid());
	g_aux.index.reserve(dim_.threads());
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for sorting buffers");
    }

    // allocate global device memory for binary mixture sampling
    for (size_t n = 0, i = 0; n < npart; n += mpart[i], ++i) {
	cuda::config dim((mpart[i] + threads() - 1) / threads(), threads());
	g_aux.index.reserve(n + dim.threads());
	dim_sample.push_back(dim);
    }

    try {
	radix_.resize(npart, dim_.threads_per_block());
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for radix sort");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::state(host_sample_type& sample, float_type box)
{
    _Base::state(sample, box, h_part.r, h_part.v);

    try {
	cuda::copy(h_part.r, g_part.r, stream_);
	cuda::copy(h_part.v, g_part.v, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const&) {
	throw exception("failed to copy phase space sample to GPU");
    }

    assign_positions();
    assign_velocities();
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::lattice()
{
    // place particles on an fcc lattice
    _Base::lattice(g_part.r);
    // randomly permute particle coordinates for binary mixture
    _Base::random_permute(g_part.r);

    assign_positions();
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::temperature(float_type temp)
{
    LOG("initialising velocities from Boltzmann distribution at temperature: " << temp);

    try {
	event_[0].record(stream_);
	_Base::boltzmann(g_part.v, temp, stream_);
	event_[1].record(stream_);
	event_[1].synchronize();
    }
    catch (cuda::error const&) {
	throw exception("failed to compute Boltzmann distributed velocities on GPU");
    }
    m_times["boltzmann"] += event_[1] - event_[0];

    assign_velocities();
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
	event_[9].record(stream_);

#ifdef USE_HILBERT_ORDER
	try {
	    permutation(stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to generate permutation for phase space sampling");
	}
#endif
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

    // heat bath coupling
    if (thermostat_steps && ++thermostat_count > thermostat_steps) {
	try {
	    _Base::boltzmann(g_part.v, thermostat_temp, stream_);
	}
	catch (cuda::error const&) {
	    throw exception("failed to compute Boltzmann distributed velocities on GPU");
	}
    }
    event_[10].record(stream_);

    // potential energy sum calculation
    try {
	reduce_en(g_part.en, stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream potential energy sum calculation on GPU");
    }
    event_[7].record(stream_);

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
    m_times["potential_energy"] += event_[7] - event_[10];
    // GPU time for maximum velocity calculation
    m_times["maximum_velocity"] += event_[0] - event_[7];

    if (thermostat_steps && thermostat_count > thermostat_steps) {
	// reset MD steps since last heatbath coupling
	thermostat_count = 0;
	// GPU time for Maxwell-Boltzmann distribution
	m_times["boltzmann"] += event_[10] - event_[6];
    }

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
	m_times["update_neighbours"] += event_[9] - event_[4];
#if defined(USE_HILBERT_ORDER)
	// GPU time for permutation sort
	m_times["permutation"] += event_[5] - event_[9];
#endif
    }

    if (!std::isfinite(reduce_en.value())) {
	throw exception("potential energy diverged");
    }

    // add to sum over maximum velocity magnitudes since last cell lists update
    v_max_sum += reduce_v_max.value();
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::sample(host_sample_type& sample) const
{
    typedef typename host_sample_type::value_type sample_type;
    typedef typename sample_type::position_sample_vector position_sample_vector;
    typedef typename sample_type::position_sample_ptr position_sample_ptr;
    typedef typename sample_type::velocity_sample_vector velocity_sample_vector;
    typedef typename sample_type::velocity_sample_ptr velocity_sample_ptr;

    static cuda::event ev0, ev1;
    static cuda::stream stream;

    try {
	ev0.record(stream);
	cuda::copy(g_part.r, h_part.r, stream);
	cuda::copy(g_part.R, h_part.R, stream);
	cuda::copy(g_part.v, h_part.v, stream);
	cuda::copy(g_part.tag, h_part.tag, stream);
	ev1.record(stream);
	ev1.synchronize();
    }
    catch (cuda::error const&) {
	throw exception("failed to copy MD simulation step results from GPU to host");
    }
    m_times["sample_memcpy"] += ev1 - ev0;

    // allocate memory for phase space sample
    for (size_t n = 0, i = 0; n < npart; n += mpart[i], ++i) {
	position_sample_ptr r(new position_sample_vector(mpart[i]));
	velocity_sample_ptr v(new velocity_sample_vector(mpart[i]));
	sample.push_back(sample_type(r, v));
    }

    // copy particle positions and velocities in binary mixture
    for (size_t i = 0, tag; i < h_part.tag.size(); ++i) {
	if ((tag = h_part.tag[i]) != gpu::VIRTUAL_PARTICLE) {
	    unsigned int const type = (tag >= mpart[0]);
	    unsigned int const n = type ? (tag - mpart[0]) : tag;
	    (*sample[type].r)[n] = h_part.r[i] + box_ * static_cast<vector_type>(h_part.R[i]);
	    (*sample[type].v)[n] = h_part.v[i];
	}
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::sample(gpu_sample_type& sample) const
{
    typedef typename gpu_sample_type::value_type sample_type;
    typedef typename sample_type::position_sample_vector position_sample_vector;
    typedef typename sample_type::position_sample_ptr position_sample_ptr;
    typedef typename sample_type::velocity_sample_vector velocity_sample_vector;
    typedef typename sample_type::velocity_sample_ptr velocity_sample_ptr;

    static cuda::event ev0, ev1;
    static cuda::stream stream;

    ev0.record(stream_);

    for (size_t n = 0, i = 0; n < npart; n += mpart[i], ++i) {
	// allocate global device memory for phase space sample
	position_sample_ptr r(new position_sample_vector(mpart[i]));
	velocity_sample_ptr v(new velocity_sample_vector(mpart[i]));
	sample.push_back(sample_type(r, v));
	// allocate additional memory to match CUDA grid dimensions
	r->reserve(dim_sample[i].threads());
	v->reserve(dim_sample[i].threads());
	// order particles by permutation
	cuda::configure(dim_sample[i].grid, dim_sample[i].block, stream);
	_gpu::sample(g_aux.index.data() + n, *r, *v);
    }

    ev1.record(stream_);
    ev1.synchronize();
    m_times["sample"] += ev1 - ev0;
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::sample(energy_sample_type& sample) const
{
    static cuda::event ev0, ev1;
    static cuda::stream stream;

    // mean potential energy per particle
    sample.en_pot = reduce_en.value() / npart;

    // mean virial equation sum per particle
    try {
	ev0.record(stream);
	reduce_virial(g_part.virial, stream);
	ev1.record(stream);
	ev1.synchronize();
	m_times["virial_sum"] += ev1 - ev0;
    }
    catch (cuda::error const& e) {
	throw exception("failed to calculate virial equation sum on GPU");
    }
    sample.virial = reduce_virial.value() / npart;

    // mean squared velocity per particle
    try {
	ev0.record(stream);
	reduce_squared_velocity(g_part.v, stream);
	ev1.record(stream);
	ev1.synchronize();
	m_times["reduce_squared_velocity"] += ev1 - ev0;
    }
    catch (cuda::error const& e) {
	throw exception("failed to calculate mean squared velocity on GPU");
    }
    sample.vv = reduce_squared_velocity.value() / npart;

    // mean velocity per particle
    try {
	ev0.record(stream);
	reduce_velocity(g_part.v, stream);
	ev1.record(stream);
	ev1.synchronize();
	m_times["reduce_velocity"] += ev1 - ev0;
    }
    catch (cuda::error const& e) {
	throw exception("failed to calculate mean velocity on GPU");
    }
    sample.v_cm = reduce_velocity.value() / npart;
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::assign_positions()
{
    // assign ascending particle numbers
    _Base::init_tags(g_part.r, g_part.tag);

    try {
	// set periodic box traversal vectors to zero
	cuda::memset(g_part.R, 0);
#ifdef USE_HILBERT_ORDER
	// order particles after Hilbert space-filling curve
	hilbert_order(stream_);
#endif
	// copy particles tags from GPU to host
	cuda::copy(g_part.tag, h_part.tag, stream_);
	// assign particles to cells
	assign_cells(stream_);
	// update neighbour lists
	update_neighbours(stream_);
	// generate permutation for phase space sampling
	permutation(stream_);
	// calculate forces
	update_forces(stream_);
	// calculate potential energy
	reduce_en(g_part.en, stream_);
	// calculate virial equation sum
	reduce_virial(g_part.virial, stream_);

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to assign particle positions on GPU");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::assign_velocities()
{
    try {
	reduce_v_max(g_part.v, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to assign particle velocities on GPU");
    }
    v_max_sum = reduce_v_max.value();
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
    _Base::update_forces(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::assign_cells(cuda::stream& stream)
{
    // compute cell indices for particle positions
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::compute_cell(g_part.r, g_aux.cell);

    // generate permutation
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::gen_index(g_aux.index);
    radix_(g_aux.cell, g_aux.index, stream);

    // compute global cell offsets in sorted particle list
    cuda::memset(g_aux.offset, 0xFF);
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::find_cell_offset(g_aux.cell, g_aux.offset);

    // assign particles to cells
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream);
    _gpu::assign_cells(g_aux.cell, g_aux.offset, g_aux.index, g_cell);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::update_neighbours(cuda::stream& stream)
{
    // mark neighbour list placeholders as virtual particles
    cuda::memset(g_nbl, 0xFF);
    // build neighbour lists
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream);
    _gpu::update_neighbours(g_cell);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::hilbert_order(cuda::stream& stream)
{
    // compute Hilbert space-filling curve for particles
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::hilbert<dimension>::curve(g_part.r, g_aux.cell);

    // generate permutation
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::gen_index(g_aux.index);
    radix_(g_aux.cell, g_aux.index, stream);

    // order particles by permutation
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::order_particles(g_aux.index, g_part_buf.r, g_part_buf.R, g_part_buf.v, g_part.tag);
    cuda::copy(g_part_buf.r, g_part.r, stream);
    cuda::copy(g_part_buf.R, g_part.R, stream);
    cuda::copy(g_part_buf.v, g_part.v, stream);
}

/**
 * generate permutation for phase space sampling
 */
template <int dimension>
void ljfluid<ljfluid_impl_gpu_neighbour<dimension> >::permutation(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, stream_);
    _gpu::gen_index(g_aux.index);
#ifdef USE_HILBERT_ORDER
    cuda::copy(g_part.tag, g_aux.cell, stream_);
    radix_(g_aux.cell, g_aux.index, stream_);
#endif
}

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
