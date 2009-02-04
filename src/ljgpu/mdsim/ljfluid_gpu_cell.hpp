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
    typedef typename _Base::host_sample_type host_sample_type;
    typedef typename _Base::gpu_sample_type gpu_sample_type;
    typedef typename _Base::energy_sample_type energy_sample_type;

    /** static implementation properties */
    typedef boost::false_type has_trajectory_gpu_sample;
    typedef boost::true_type has_energy_gpu_sample;

public:
    /** set number of CUDA execution threads */
    void threads(unsigned int value);
    /** set desired average cell occupancy */
    void cell_occupancy(float_type value);

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
    /** sample thermodynamic equilibrium properties */
    void sample(energy_sample_type& sample) const;

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
    boost::array<cuda::event, 7> event_;

    using _Base::reduce_squared_velocity;
    using _Base::reduce_velocity;
    using _Base::reduce_en;
    using _Base::reduce_virial;

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
    /** maximum absolute velocity */
    reduce<tag::max, float> reduce_v_max;
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
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::state(host_sample_type& sample, float_type box)
{
    cuda::host::vector<float4> h_r(npart);
    cuda::host::vector<gpu_vector_type> h_v(npart);
    _Base::state(sample, box, h_r, h_v);

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
    LOG("initialising velocities from Boltzmann distribution at temperature: " << temp);

    cuda::vector<gpu_vector_type> g_v(npart);
    cuda::host::vector<gpu_vector_type> h_v(npart);
    g_v.reserve(dim_.threads());
    h_v.reserve(dim_.threads());

    try {
	event_[0].record(stream_);
	_Base::boltzmann(g_v, temp, stream_);
	event_[1].record(stream_);
	event_[1].synchronize();
    }
    catch (cuda::error const&) {
	throw exception("failed to compute Boltzmann distributed velocities on GPU");
    }
    m_times["boltzmann"] += event_[1] - event_[0];

    cuda::copy(g_v, h_v);
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

    // maximum velocity calculation
    try {
	reduce_v_max(g_part.v, stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream maximum velocity calculation on GPU");
    }
    event_[3].record(stream_);
    event_[3].synchronize();

    v_max_sum += reduce_v_max.value();

    // update cell lists
    if (v_max_sum * timestep_ > r_skin / 2) {
	try {
	    update_cells(stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to stream cell list update on GPU");
	}
	event_[4].record(stream_);

	try {
	    copy_cells(stream_);
	}
	catch (cuda::error const& e) {
	    throw exception("failed to replicate cell lists on GPU");
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

    m_times["mdstep"] += event_[0] - event_[1];
    m_times["velocity_verlet"] += event_[2] - event_[1];
    m_times["maximum_velocity"] += event_[3] - event_[2];

    if (v_max_sum * timestep_ > r_skin / 2) {
	// reset sum over maximum velocity magnitudes to zero
	v_max_sum = 0;

	m_times["update_cells"] += event_[4] - event_[3];
	m_times["memcpy_cells"] += event_[5] - event_[4];
    }

    m_times["update_forces"] += event_[6] - event_[5];
    m_times["potential_energy"] += event_[0] - event_[6];

    if (!std::isfinite(reduce_en.value())) {
	throw potential_energy_divergence();
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::sample(host_sample_type& sample) const
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
    catch (cuda::error const& e) {
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
    size_t count = 0;
    for (size_t i = 0, tag; i < h_part.tag.size(); ++i) {
	if ((tag = h_part.tag[i]) != gpu::VIRTUAL_PARTICLE) {
	    unsigned int const type = (tag >= mpart[0]);
	    unsigned int const n = type ? (tag - mpart[0]) : tag;
	    (*sample[type].r)[n] = h_part.r[i] + box_ * static_cast<vector_type>(h_part.R[i]);
	    (*sample[type].v)[n] = h_part.v[i];
	    count++;
	}
    }
    if (count != npart) {
	throw exception("particle loss while updating cell lists");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_cell<dimension> >::sample(energy_sample_type& sample) const
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
	// potential energy sum calculation
	reduce_en(g_part.en, stream_);

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
    typename cuda::host::vector<unsigned int>::const_iterator tag;
    typename cuda::host::vector<gpu_vector_type>::iterator v;

    // assign particle velocities to cell placeholders
    for (tag = h_part.tag.begin(), v = h_part.v.begin(); tag != h_part.tag.end(); ++tag, ++v) {
	if (*tag != gpu::VIRTUAL_PARTICLE) {
	    *v = h_v[*tag];
	}
    }

    try {
	cuda::copy(h_part.v, g_part.v, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to assign velocities to cell placeholders");
    }

    v_max_sum = 0;
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
