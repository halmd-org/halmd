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
    typedef typename _Base::sample_type sample_type;
    typedef typename sample_type::sample_visitor sample_visitor;

public:
    /** set number of particles in system */
    void particles(unsigned int value);
    /** set number of CUDA execution threads */
    void threads(unsigned int value);

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

private:
    /** first leapfrog step of integration of differential equations of motion */
    void velocity_verlet(cuda::stream& stream);
    /** Lennard-Jones force calculation */
    void update_forces(cuda::stream& stream);

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

    /** CUDA events for kernel timing */
    boost::array<cuda::event, 9> event_;

    /** potential energy sum */
    reduce<tag::sum, dfloat> reduce_en;
    /** virial equation sum */
    reduce<tag::sum, dfloat> reduce_virial;

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
	/** blockwise maximum particles velocity magnitudes */
	cuda::host::vector<float> v_max;
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
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
    } g_part;

};

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::particles(unsigned int value)
{
    _Base::particles(value);

    // allocate global device memory for system state
    try {
	g_part.r.resize(npart);
	g_part.R.resize(npart);
	g_part.v.resize(npart);
	g_part.f.resize(npart);
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
	g_part.en.reserve(dim_.threads());
	g_part.virial.reserve(dim_.threads());
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate global device memory for placeholder particles");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::sample(sample_visitor visitor)
{
    _Base::sample(visitor);

    try {
	// copy periodically reduced particle positions from host to GPU
	std::copy(m_sample[0].r.begin(), m_sample[0].r.end(), h_part.r.begin());
	cuda::copy(h_part.r, g_part.r, stream_);
	// set periodic box traversal vectors to zero
	cuda::memset(g_part.R, 0);
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

	// wait for GPU operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to restore system state from phase space sample");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::lattice()
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
	// compute particle lattice positions on GPU
	event_[0].record(stream_);
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::lattice::fcc(g_part.r, n, box_);
	event_[1].record(stream_);
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
void ljfluid<ljfluid_impl_gpu_square<dimension> >::temperature(float_type temp)
{
    LOG("initialising velocities from Maxwell-Boltzmann distribution at temperature: " << temp);
    boltzmann(g_part.v, h_part.v, temp);
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

    if (!std::isfinite((double) reduce_en.value())) {
	throw exception("potential energy diverged");
    }
}

template <int dimension>
void ljfluid<ljfluid_impl_gpu_square<dimension> >::copy()
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
	event_[0].record(stream_);

	// wait for CUDA operations to finish
	event_[0].synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy MD simulation step results from GPU to host");
    }

    for (unsigned int j = 0; j < npart; ++j) {
	// copy periodically extended particle positions
	m_sample[0].r[j] = h_part.r[j] + (vector_type) h_part.R[j] * box_;
	// copy particle velocities
	m_sample[0].v[j] = h_part.v[j];
    }

    // mean potential energy per particle
    m_sample.en_pot = (double) reduce_en.value() / npart;
    // mean virial equation sum per particle
    m_sample.virial = (double) reduce_virial.value() / npart;

    // GPU time for sample memcpy
    m_times["sample_memcpy"] += event_[0] - event_[1];
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
    cuda::configure(dim_.grid, dim_.block, dim_.threads_per_block() * sizeof(gpu_vector_type), stream);
    if (r_smooth > 0 && thermostat_nu > 0) {
	_gpu::mdstep_smooth_nvt(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
    }
    else if (r_smooth > 0) {
	_gpu::mdstep_smooth(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
    }
    else if (thermostat_nu > 0) {
	_gpu::mdstep_nvt(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
    }
    else {
	_gpu::mdstep(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
    }
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_GPU_SQUARE_HPP */
