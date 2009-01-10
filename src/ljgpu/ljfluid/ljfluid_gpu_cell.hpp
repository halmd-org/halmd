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

#ifndef LJGPU_LJFLUID_LJFLUID_GPU_CELL_HPP
#define LJGPU_LJFLUID_LJFLUID_GPU_CELL_HPP

#include <ljgpu/ljfluid/base_gpu.hpp>
#include <ljgpu/ljfluid/gpu/lattice.hpp>
#include <ljgpu/math/stat.hpp>

namespace ljgpu
{

template <typename ljfluid_impl>
class ljfluid;

template <int dimension>
class ljfluid<ljfluid_impl_gpu_cell<dimension> >
    : public ljfluid_base_gpu<ljfluid_impl_gpu_cell<dimension> >
{
public:
    typedef ljfluid_base_gpu<ljfluid_impl_gpu_cell<dimension> > _Base;
    typedef gpu::ljfluid<ljfluid_impl_gpu_cell<dimension> > _gpu;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::gpu_vector_type gpu_vector_type;
    typedef typename _Base::sample_type sample_type;
    typedef typename sample_type::sample_visitor sample_visitor;

public:
    /** initialise fluid from program options */
    ljfluid(options const& opt);

    using _Base::particles;
    using _Base::density;
    using _Base::box;
    using _Base::timestep;
    using _Base::cutoff_radius;
#ifdef USE_POTENTIAL_SMOOTHING
    using _Base::potential_smoothing;
#endif

    /** set number of CUDA execution threads */
    void threads(unsigned int value);
    /** set desired average cell occupancy */
    void cell_occupancy(float_type value);

    /** restore system state from phase space sample */
    void restore(sample_visitor visitor);
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

    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;

private:
    /** first leapfrog step of integration of differential equations of motion */
    void velocity_verlet(cuda::stream& stream);
    /** Lennard-Jones force calculation */
    void update_forces(cuda::stream& stream);
    /** update cell lists */
    void update_cells(cuda::stream& stream);
    /** replicate cell lists */
    void copy_cells(cuda::stream& stream);

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

    using _Base::dim_;
    using _Base::stream_;
    using _Base::rng_;

    /** CUDA execution dimensions for cell-specific kernels */
    cuda::config dim_cell_;
    /** CUDA events for kernel timing */
    boost::array<cuda::event, 9> event_;

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
	/** blockwise maximum particles velocity magnitudes */
	cuda::host::vector<float> v_max;
	/** potential energies per particle */
	cuda::host::vector<float> en;
	/** virial equation sums per particle */
	cuda::host::vector<float> virial;
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

    /** system state double buffer in global device memory */
    struct {
	/** periodically reduced particle positions */
	cuda::vector<gpu_vector_type> r;
	/** periodically extended particle positions */
	cuda::vector<gpu_vector_type> R;
	/** particle velocities */
	cuda::vector<gpu_vector_type> v;
	/** particle tags */
	cuda::vector<int> tag;
    } g_part_buf;
};

template <int dimension>
ljfluid<ljfluid_impl_gpu_cell<dimension> >::ljfluid(options const& opt)
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
    cell_occupancy(opt["cell-occupancy"].as<float>());
    threads(opt["threads"].as<unsigned int>());
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::cell_occupancy(float_type value)
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

    try {
	cuda::copy(ncell, _gpu::ncell);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy cell parameters to device symbols");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::threads(unsigned int value)
{
    _Base::threads(value);

    // set CUDA execution dimensions for cell-specific kernels
    dim_cell_ = cuda::config(dim3(powf(ncell, dimension)), dim3(cell_size_));
    LOG("number of cell CUDA execution blocks: " << dim_cell_.blocks_per_grid());
    LOG("number of cell CUDA execution threads: " << dim_cell_.threads_per_block());

    // allocate page-locked host memory for cell placeholders
    try {
	h_part.r.resize(dim_cell_.threads());
	h_part.R.resize(dim_cell_.threads());
	h_part.v.resize(dim_cell_.threads());
	h_part.tag.resize(dim_cell_.threads());
	// particle forces reside only in GPU memory
	h_part.en.resize(dim_cell_.threads());
	h_part.virial.resize(dim_cell_.threads());
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate page-locked host memory cell placeholders");
    }

    // allocate global device memory for cell placeholders
    try {
	g_part.r.resize(dim_cell_.threads());
	g_part.R.resize(dim_cell_.threads());
	g_part.v.resize(dim_cell_.threads());
	g_part.tag.resize(dim_cell_.threads());
	g_part.f.resize(dim_cell_.threads());
	g_part.en.resize(dim_cell_.threads());
	g_part.virial.resize(dim_cell_.threads());
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory cell placeholders");
    }

    // allocate global device memory for double buffers
    try {
	g_part_buf.r.resize(dim_cell_.threads());
	g_part_buf.R.resize(dim_cell_.threads());
	g_part_buf.v.resize(dim_cell_.threads());
	g_part_buf.tag.resize(dim_cell_.threads());
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory double buffers");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::restore(sample_visitor visitor)
{
    // read phase space sample
    visitor(m_sample.r, m_sample.v);

    try {
	// copy periodically reduced particle positions from host to GPU
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.r[i] = make_float(m_sample.r[i]);
	}
	cuda::copy(h_part.r, g_part.R, stream_);
	// assign particles to cells
	event_[0].record(stream_);
	cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	_gpu::assign_cells(g_part.R, g_part.r, g_part.tag);
	event_[1].record(stream_);
	// replicate to periodically extended particle positions
	cuda::copy(g_part.r, g_part.R, stream_);
	// calculate forces
	update_forces(stream_);

	// copy particle number tags from GPU to host
	cuda::copy(g_part.tag, h_part.tag, stream_);
	// wait for CUDA operations to finish
	stream_.synchronize();

	// assign velocities to cell placeholders
	float vv_max = 0;
	for (unsigned int i = 0; i < nplace; ++i) {
	    // particle number
	    const int tag = h_part.tag[i];
	    if (tag != _gpu::VIRTUAL_PARTICLE) {
		h_part.v[i] = make_float(m_sample.v[tag]);
		// calculate maximum squared velocity
		vv_max = std::max(vv_max, m_sample.v[tag] * m_sample.v[tag]);
	    }
	}
	// set initial sum over maximum velocity magnitudes since last cell lists update
	v_max_sum = std::sqrt(vv_max);
	cuda::copy(h_part.v, g_part.v, stream_);

	// wait for GPU operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to restore system state from phase space sample");
    }

    m_times["init_cells"] += event_[1] - event_[0];
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::lattice()
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
	g_part.r.reserve(dim_.threads());
	// compute particle lattice positions on GPU
	event_[0].record(stream_);
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::lattice::fcc(g_part.r, n, box_);
	event_[1].record(stream_);
	// TODO randomly permute particles to increase force summing accuracy
	cuda::copy(g_part.r, g_part_buf.r, stream_);
	// assign particles to cells
	event_[2].record(stream_);
	cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	_gpu::assign_cells(g_part_buf.r, g_part.r, g_part.tag);
	event_[3].record(stream_);
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0;
	// replicate particle positions to periodically extended positions
	cuda::copy(g_part.r, g_part.R, stream_);
	// calculate forces, potential energy and virial equation sum
	cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	_gpu::mdstep(g_part.r, g_part.v, g_part.f, g_part.tag, g_part.en, g_part.virial);

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute particle lattice positions on GPU");
    }

    // CUDA time for lattice generation
    m_times["lattice"] += event_[1] - event_[0];
    // CUDA time for cell lists initialisation
    m_times["init_cells"] += event_[3] - event_[2];
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::temperature(float_type temp)
{
    LOG("initialising velocities from Maxwell-Boltzmann distribution at temperature: " << temp);
    try {
	g_part.v.reserve(dim_.threads());
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
	// copy particle number tags from GPU to host
	cuda::copy(g_part.tag, h_part.tag, stream_);
	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute Maxwell-Boltzmann distributed velocities on GPU");
    }

    // CUDA time for Maxwell-Boltzmann distribution
    m_times["boltzmann"] += event_[1] - event_[0];

    // compute center of mass velocity
    vector_type v_cm = mean(m_sample.v.begin(), m_sample.v.end());
    // set center of mass velocity to zero
    for (unsigned int i = 0; i < npart; ++i) {
	m_sample.v[i] -= v_cm;
    }

    try {
	// maximum squared velocity
	float vv_max = 0;
	// assign velocities to cell placeholders
	for (unsigned int i = 0; i < nplace; ++i) {
	    // particle number
	    const int n = h_part.tag[i];
	    if (n != _gpu::VIRTUAL_PARTICLE) {
		// assign velocity to cell placeholder
		h_part.v[i] = make_float(m_sample.v[n]);
		// calculate maximum squared velocity
		vv_max = std::max(vv_max, m_sample.v[n] * m_sample.v[n]);
	    }
	}
	// initialise sum over maximum velocity magnitudes since last cell lists update
	v_max_sum = std::sqrt(vv_max);
	// copy particle velocities from host to GPU
	cuda::copy(h_part.v, g_part.v, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to set center of mass velocity to zero");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::stream()
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
	try {
	    update_cells(stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to stream cell list update on GPU");
	}
	event_[3].record(stream_);

	try {
	    copy_cells(stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to replicate cell lists on GPU");
	}
    }
    event_[4].record(stream_);

    // Lennard-Jones force calculation
    try {
	update_forces(stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream force calculation on GPU");
    }
    event_[0].record(stream_);
}

/**
 * synchronize MD simulation step on GPU
 */
template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::mdstep()
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
    // CUDA time for Lennard-Jones force update
    m_times["update_forces"] += event_[0] - event_[4];

    if (v_max_sum * timestep_ > r_skin / 2) {
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0;
	// CUDA time for cell lists update
	m_times["update_cells"] += event_[3] - event_[2];
	// CUDA time for cell lists memcpy
	m_times["memcpy_cells"] += event_[4] - event_[3];
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::copy()
{
    // mean potential energy per particle
    m_sample.en_pot = 0;
    // mean virial equation sum per particle
    m_sample.virial = 0;

    // copy MD simulation step results from GPU to host
    try {
	event_[1].record(stream_);
	// copy periodically reduce particles positions
	cuda::copy(g_part.r, h_part.r, stream_);
	// copy periodically extended particles positions
	cuda::copy(g_part.R, h_part.R, stream_);
	// copy particle velocities
	cuda::copy(g_part.v, h_part.v, stream_);
	// copy particle number tags
	cuda::copy(g_part.tag, h_part.tag, stream_);
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

    // number of particles found in cells
    unsigned int count = 0;
    // maximum squared velocity
    float vv_max = 0;

    for (unsigned int i = 0; i < nplace; ++i) {
	// particle number
	const int n = h_part.tag[i];
	// check if real particle
	if (n != _gpu::VIRTUAL_PARTICLE) {
	    // copy periodically reduced particle positions
	    m_sample.r[n] = h_part.r[i];
	    // copy periodically extended particle positions
	    m_sample.R[n] = h_part.R[i];
	    // copy particle velocities
	    m_sample.v[n] = h_part.v[i];
	    // calculate mean potential energy per particle
	    m_sample.en_pot += (h_part.en[i] - m_sample.en_pot) / ++count;
	    // calculate mean virial equation sum per particle
	    m_sample.virial += (h_part.virial[i] - m_sample.virial) / count;
	    // calculate maximum squared velocity
	    vv_max = std::max(vv_max, m_sample.v[n] * m_sample.v[n]);
	}
    }
    // add to sum over maximum velocity magnitudes since last cell lists update
    v_max_sum += std::sqrt(vv_max);
    // validate number of particles
    if (count != npart) {
	throw exception("particle loss while updating cell lists");
    }

    // ensure that system is still in valid state after MD step
    if (std::isnan(m_sample.en_pot)) {
	throw exception("potential energy diverged due to excessive timestep or density");
    }

    // GPU time for sample memcpy
    m_times["sample_memcpy"] += event_[0] - event_[1];
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::velocity_verlet(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::inteq(g_part.r, g_part.R, g_part.v, g_part.f);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::update_forces(cuda::stream& stream)
{
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
    _gpu::mdstep(g_part.r, g_part.v, g_part.f, g_part.tag, g_part.en, g_part.virial);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::update_cells(cuda::stream& stream)
{
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream);
    _gpu::update_cells(g_part.r, g_part.R, g_part.v, g_part.tag, g_part_buf.r, g_part_buf.R, g_part_buf.v, g_part_buf.tag);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::copy_cells(cuda::stream& stream)
{
    cuda::copy(g_part_buf.r, g_part.r, stream);
    cuda::copy(g_part_buf.R, g_part.R, stream);
    cuda::copy(g_part_buf.v, g_part.v, stream);
    cuda::copy(g_part_buf.tag, g_part.tag, stream);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::attrs(H5::Group const& param) const
{
    _Base::attrs(param);

    H5xx::group node(param.openGroup("mdsim"));
    node["cells"] = cells();
    node["placeholders"] = placeholders();
    node["cell_length"] = cell_length();
    node["cell_occupancy"] = cell_occupancy();
}

} // namespace ljgpu

#endif /* ! LJGPU_LJFLUID_LJFLUID_GPU_CELL_HPP */
