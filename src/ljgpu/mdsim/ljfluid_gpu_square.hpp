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

#ifndef LJGPU_MDSIM_LJFLUID_GPU_SQUARE_HPP
#define LJGPU_MDSIM_LJFLUID_GPU_SQUARE_HPP

#include <ljgpu/algorithm/reduce.hpp>
#include <ljgpu/mdsim/ljfluid_gpu_base.hpp>
#include <ljgpu/mdsim/gpu/lattice.hpp>

namespace ljgpu
{

template <typename ljfluid_impl>
class ljfluid;

template <int dimension>
class ljfluid<ljfluid_impl_gpu_square<dimension> >
    : public ljfluid_gpu_base<ljfluid_impl_gpu_square<dimension> >
{
public:
    typedef ljfluid_gpu_base<ljfluid_impl_gpu_square<dimension> > _Base;
    typedef gpu::ljfluid<ljfluid_impl_gpu_square<dimension> > _gpu;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::gpu_vector_type gpu_vector_type;
    typedef typename _Base::host_sample_type host_sample_type;
    typedef typename _Base::gpu_sample_type gpu_sample_type;
    typedef typename _Base::energy_sample_type energy_sample_type;

    /** static implementation properties */
    typedef boost::true_type has_trajectory_gpu_sample;
    typedef boost::true_type has_energy_gpu_sample;

public:
    /** set number of particles in system */
    template <typename T>
    void particles(T const& value);
    /** set number of CUDA execution threads */
    void threads(unsigned int value);

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

private:
    /** assign particle positions */
    void assign_positions();
    /** first leapfrog step of integration of differential equations of motion */
    void velocity_verlet(cuda::stream& stream);
    /** Lennard-Jones force calculation */
    void update_forces(cuda::stream& stream);

private:
    using _Base::box_;
    using _Base::density_;
    using _Base::dim_;
    using _Base::ensemble_;
    using _Base::m_times;
    using _Base::mixture_;
    using _Base::mpart;
    using _Base::npart;
    using _Base::potential_;
    using _Base::r_cut;
    using _Base::stream_;
    using _Base::timestep_;

    /** CUDA events for kernel timing */
    boost::array<cuda::event, 9> event_;
    /** CUDA execution dimensions for phase space sampling */
    std::vector<cuda::config> dim_sample;

    using _Base::reduce_squared_velocity;
    using _Base::reduce_velocity;
    using _Base::reduce_en;
    using _Base::reduce_virial;
    using _Base::reduce_v_max;

    /** system state in page-locked host memory */
    struct {
	/** tagged periodically reduced particle positions */
	cuda::host::vector<float4> mutable r;
	/** periodic box traversal vectors */
	cuda::host::vector<gpu_vector_type> mutable R;
	/** particle velocities */
	cuda::host::vector<gpu_vector_type> mutable v;
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

};

template <int dimension>
template <typename T>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::particles(T const& value)
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
    // allocate page-locked host memory for system state
    try {
	h_part.r.resize(npart);
	h_part.R.resize(npart);
	h_part.v.resize(npart);
	// particle forces reside only in GPU memory
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate page-locked host memory for system state");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::threads(unsigned int value)
{
    _Base::threads(value);

    // allocate global device memory for placeholder particles
    try {
	g_part.r.reserve(dim_.threads());
	g_part.R.reserve(dim_.threads());
	g_part.v.reserve(dim_.threads());
	g_part.f.reserve(dim_.threads());
	g_part.tag.reserve(dim_.threads());
	g_part.en.reserve(dim_.threads());
	g_part.virial.reserve(dim_.threads());
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for placeholder particles");
    }

    // allocate global device memory for binary mixture sampling
    for (size_t n = 0, i = 0; n < npart; n += mpart[i], ++i) {
	cuda::config dim((mpart[i] + threads() - 1) / threads(), threads());
	g_part.r.reserve(n + dim.threads());
	g_part.R.reserve(n + dim.threads());
	g_part.v.reserve(n + dim.threads());
	dim_sample.push_back(dim);
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::state(host_sample_type& sample, float_type box)
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
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::lattice()
{
    // place particles on an fcc lattice
    _Base::lattice(g_part.r);
    // randomly permute particle coordinates for binary mixture
    _Base::random_permute(g_part.r);

    assign_positions();
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::temperature(float_type temp)
{
    _Base::boltzmann(g_part.v, h_part.v, temp);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::stream()
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
	reduce_en(g_part.en, stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream potential energy sum calculation on GPU");
    }
    event_[4].record(stream_);

    // virial equation sum calculation
    try {
	reduce_virial(g_part.virial, stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream virial equation sum calculation on GPU");
    }
    event_[0].record(stream_);
}

/**
 * synchronize MD simulation step on GPU
 */
template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::mdstep()
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
    m_times["update_forces"] += event_[3] - event_[2];
    // GPU time for potential energy sum calculation
    m_times["potential_energy"] += event_[4] - event_[3];
    // GPU time for virial equation sum calculation
    m_times["virial_sum"] += event_[0] - event_[4];

    if (!std::isfinite(reduce_en.value())) {
	throw exception("potential energy diverged");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::sample(host_sample_type& sample) const
{
    typedef typename host_sample_type::position_sample_vector position_sample_vector;
    typedef typename host_sample_type::position_sample_ptr position_sample_ptr;
    typedef typename host_sample_type::velocity_sample_vector velocity_sample_vector;
    typedef typename host_sample_type::velocity_sample_ptr velocity_sample_ptr;

    static cuda::event ev0, ev1;
    static cuda::stream stream;

    try {
	ev0.record(stream);
	cuda::copy(g_part.r, h_part.r, stream);
	cuda::copy(g_part.R, h_part.R, stream);
	cuda::copy(g_part.v, h_part.v, stream);
	ev1.record(stream);
	ev1.synchronize();
    }
    catch (cuda::error const&) {
	throw exception("failed to copy MD simulation step results from GPU to host");
    }
    m_times["sample_memcpy"] += ev1 - ev0;

    for (size_t n = 0, i = 0; n < npart; ++i) {
	// allocate memory for trajectory sample
	position_sample_ptr r(new position_sample_vector);
	velocity_sample_ptr v(new velocity_sample_vector);
	sample.r.push_back(r);
	sample.v.push_back(v);
	r->reserve(mpart[i]);
	v->reserve(mpart[i]);
	// assign particle positions and velocities of homogenous type
	for (size_t j = 0; j < mpart[i]; ++j, ++n) {
	    r->push_back(h_part.r[n] + box_ * static_cast<vector_type>(h_part.R[n]));
	    v->push_back(h_part.v[n]);
	}
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::sample(gpu_sample_type& sample) const
{
    typedef typename gpu_sample_type::position_sample_vector position_sample_vector;
    typedef typename gpu_sample_type::position_sample_ptr position_sample_ptr;
    typedef typename gpu_sample_type::velocity_sample_vector velocity_sample_vector;
    typedef typename gpu_sample_type::velocity_sample_ptr velocity_sample_ptr;

    static cuda::event ev0, ev1;
    static cuda::stream stream;

    ev0.record(stream_);

    for (size_t n = 0, i = 0; n < npart; n += mpart[i], ++i) {
	// allocate global device memory for phase space sample
	position_sample_ptr r(new position_sample_vector(mpart[i]));
	velocity_sample_ptr v(new velocity_sample_vector(mpart[i]));
	sample.r.push_back(r);
	sample.v.push_back(v);
	// allocate additional memory to match CUDA grid dimensions
	r->reserve(dim_sample[i].threads());
	v->reserve(dim_sample[i].threads());
	// sample trajectories
	cuda::configure(dim_sample[i].grid, dim_sample[i].block, stream);
	_gpu::sample(g_part.r.data() + n, g_part.R.data() + n, g_part.v.data() + n, *r, *v);
    }

    ev1.record(stream_);
    ev1.synchronize();
    m_times["sample"] += ev1 - ev0;
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::sample(energy_sample_type& sample) const
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
void ljfluid<ljfluid_impl_gpu_square<dimension> >::assign_positions()
{
    // assign ascending particle numbers
    _Base::init_tags(g_part.r, g_part.tag);

    try {
	// set periodic box traversal vectors to zero
	cuda::memset(g_part.R, 0);
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
void ljfluid<ljfluid_impl_gpu_square<dimension> >::velocity_verlet(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, stream);
    _gpu::inteq(g_part.r, g_part.R, g_part.v, g_part.f);
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::update_forces(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, dim_.threads_per_block() * (dimension + 1) * sizeof(int), stream);
    _Base::update_forces(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_GPU_SQUARE_HPP */
