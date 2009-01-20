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

#ifndef LJGPU_MDSIM_LJFLUID_GPU_CELL_HPP
#define LJGPU_MDSIM_LJFLUID_GPU_CELL_HPP

#include <ljgpu/mdsim/ljfluid_gpu_base.hpp>
#include <ljgpu/mdsim/gpu/lattice.hpp>
#include <ljgpu/math/stat.hpp>

namespace ljgpu
{

template <typename ljfluid_impl>
class ljfluid;

template <int dimension>
class ljfluid<ljfluid_impl_gpu_cell<dimension> >
    : public ljfluid_gpu_base<ljfluid_impl_gpu_cell<dimension> >
{
public:
    typedef ljfluid_gpu_base<ljfluid_impl_gpu_cell<dimension> > _Base;
    typedef gpu::ljfluid<ljfluid_impl_gpu_cell<dimension> > _gpu;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::gpu_vector_type gpu_vector_type;
    typedef typename _Base::sample_type sample_type;
    typedef typename sample_type::sample_visitor sample_visitor;

public:
    /** set number of CUDA execution threads */
    void threads(unsigned int value);
    /** set desired average cell occupancy */
    void cell_occupancy(float_type value);

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

    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

private:
    /** assign particle positions */
    void assign_positions(cuda::vector<float4>& g_r);
    /** assign particle velocities */
    void assign_velocities(cuda::host::vector<gpu_vector_type>& h_v);
    /** first leapfrog step of integration of differential equations of motion */
    void velocity_verlet(cuda::stream& stream);
    /** Lennard-Jones force calculation */
    void update_forces(cuda::stream& stream);
    /** update cell lists */
    void update_cells(cuda::stream& stream);
    /** replicate cell lists */
    void copy_cells(cuda::stream& stream);

private:
    using _Base::box_;
    using _Base::density_;
    using _Base::dim_;
    using _Base::ensemble_;
    using _Base::m_sample;
    using _Base::m_times;
    using _Base::mixture_;
    using _Base::mpart;
    using _Base::npart;
    using _Base::potential_;
    using _Base::r_cut;
    using _Base::stream_;
    using _Base::timestep_;

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
	/** tagged periodically reduced particle positions */
	cuda::host::vector<float4> r;
	/** periodic box traversal vectors */
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
	/** tagged periodically reduced particle positions */
	cuda::vector<float4> r;
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

    /** system state double buffer in global device memory */
    struct {
	/** tagged periodically reduced particle positions */
	cuda::vector<float4> r;
	/** periodic box traversal vectors */
	cuda::vector<gpu_vector_type> R;
	/** particle velocities */
	cuda::vector<gpu_vector_type> v;
    } g_part_buf;
};

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::cell_occupancy(float_type value)
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
    // set cell skin
    r_skin = std::max(0.f, cell_length_ - r_cut_max);
    LOG("cell skin: " << r_skin);

    // set total number of cell placeholders
    nplace = std::pow(ncell, dimension) * cell_size_;
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
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory double buffers");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::sample(sample_visitor visitor)
{
    cuda::host::vector<float4> h_r(npart);
    cuda::host::vector<gpu_vector_type> h_v(npart);
    _Base::sample(visitor, h_r, h_v);

    cuda::vector<float4> g_r(npart);
    g_r.reserve(dim_.threads());

    try {
	cuda::copy(h_r, g_r, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const&) {
	throw exception("failed to copy particle positions to GPU");
    }

    assign_positions(g_r);
    assign_velocities(h_v);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::lattice()
{
    // place particles on an fcc lattice
    cuda::vector<float4> g_r(npart);
    g_r.reserve(dim_.threads());
    _Base::lattice(g_r);
    // randomly permute particle coordinates for binary mixture
    _Base::random_permute(g_r);

    assign_positions(g_r);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::temperature(float_type temp)
{
    cuda::vector<gpu_vector_type> g_v(npart);
    cuda::host::vector<gpu_vector_type> h_v(npart);
    g_v.reserve(dim_.threads());
    h_v.reserve(dim_.threads());
    _Base::boltzmann(g_v, h_v, temp);

    assign_velocities(h_v);
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
	// copy periodic box traversal vectors
	cuda::copy(g_part.R, h_part.R, stream_);
	// copy particle velocities
	cuda::copy(g_part.v, h_part.v, stream_);
	if (v_max_sum * timestep_ > r_skin / 2) {
	    // reset sum over maximum velocity magnitudes to zero
	    v_max_sum = 0;
	    // copy particle number tags
	    cuda::copy(g_part.tag, h_part.tag, stream_);
	}
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
	// particle tag
	int const tag = h_part.tag[i];
	// skip virtual particles
	if (tag == gpu::VIRTUAL_PARTICLE) continue;

	// A or B particle type
	unsigned int const type = (static_cast<unsigned int>(tag) >= mpart[0]);
	// particle number
	unsigned int const n = type ? (tag - mpart[0]) : tag;
	// copy periodically extended particle positions
	m_sample[type].r[n] = h_part.r[i] + (vector_type) h_part.R[i] * box_;
	// copy particle velocities
	m_sample[type].v[n] = h_part.v[i];
	// calculate mean potential energy per particle
	m_sample.en_pot += (h_part.en[i] - m_sample.en_pot) / ++count;
	// calculate mean virial equation sum per particle
	m_sample.virial += (h_part.virial[i] - m_sample.virial) / count;
	// calculate maximum squared velocity
	vv_max = std::max(vv_max, m_sample[type].v[n] * m_sample[type].v[n]);
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
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::assign_positions(cuda::vector<float4>& g_r)
{
    // assign ascending particle numbers
    _Base::init_tags(g_r, g_part.tag);

    try {
	// assign particles to cells
	event_[0].record(stream_);
	cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
	_gpu::assign_cells(g_r, g_part.r, g_part.tag);
	event_[1].record(stream_);
	// set periodic box traversal vectors to zero
	cuda::memset(g_part.R, 0);
	// copy particles tags from GPU to host
	cuda::copy(g_part.tag, h_part.tag, stream_);
	// calculate forces, potential energy and virial equation sum
	update_forces(stream_);

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to assign positions to cell placeholders");
    }

    // CUDA time for cell lists initialisation
    m_times["init_cells"] += event_[1] - event_[0];
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::assign_velocities(cuda::host::vector<gpu_vector_type>& h_v)
{
    float vv_max = 0;
    for (unsigned int i = 0; i < nplace; ++i) {
	int const tag = h_part.tag[i];
	if (tag != gpu::VIRTUAL_PARTICLE) {
	    h_part.v[i] = h_v[tag];
	    vv_max = std::max(vv_max, (vector_type) h_v[tag] * h_v[tag]);
	}
    }
    v_max_sum = std::sqrt(vv_max);

    try {
	cuda::copy(h_part.v, g_part.v, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to assign velocities to cell placeholders");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::velocity_verlet(cuda::stream& stream)
{
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream);
    _gpu::inteq(g_part.r, g_part.R, g_part.v, g_part.f);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::update_forces(cuda::stream& stream)
{
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream_);
    _Base::update_forces(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::update_cells(cuda::stream& stream)
{
    cuda::configure(dim_cell_.grid, dim_cell_.block, stream);
    _gpu::update_cells(g_part.r, g_part.R, g_part.v, g_part_buf.r, g_part_buf.R, g_part_buf.v, g_part.tag);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::copy_cells(cuda::stream& stream)
{
    cuda::copy(g_part_buf.r, g_part.r, stream);
    cuda::copy(g_part_buf.R, g_part.R, stream);
    cuda::copy(g_part_buf.v, g_part.v, stream);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::param(H5param& param) const
{
    _Base::param(param);

    H5xx::group node(param["mdsim"]);
    node["cells"] = ncell;
    node["placeholders"] = nplace;
    node["cell_length"] = cell_length_;
    node["cell_occupancy"] = cell_occupancy_;
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_GPU_CELL_HPP */
