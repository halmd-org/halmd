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

#include <boost/foreach.hpp>
#include <cuda_wrapper.hpp>
#include <ljgpu/algorithm/reduce.hpp>
#include <ljgpu/math/stat.hpp>
#include <ljgpu/mdsim/gpu/lattice.hpp>
#include <ljgpu/mdsim/gpu/ljfluid_square.hpp>
#include <ljgpu/mdsim/ljfluid_traits.hpp>
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
class ljfluid_gpu_impl_square
{
public:
    typedef typename ljfluid_gpu_traits<dimension>::float_type float_type;
    typedef typename ljfluid_gpu_traits<dimension>::vector_type vector_type;
    typedef typename ljfluid_gpu_traits<dimension>::gpu_vector_type gpu_vector_type;
    typedef typename ljfluid_gpu_traits<dimension>::trajectory_sample trajectory_sample;
    typedef typename trajectory_sample::visitor trajectory_visitor;

public:
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
    boost::array<cuda::event, 5> event_;

    /** trajectory sample in swappable host memory */
    trajectory_sample m_sample;
    /** GPU time accumulators */
    perf_counters m_times;

private:
    /** potential energy sum */
    reduce<tag::sum, dfloat> reduce_en;
    /** virial equation sum */
    reduce<tag::sum, dfloat> reduce_virial;

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
	/** potential energies per particle */
	cuda::vector<float> en;
	/** virial equation sums per particle */
	cuda::vector<float> virial;
    } g_part;

};

template <int dimension>
void ljfluid_gpu_impl_square<dimension>::particles(unsigned int value)
{
    // copy particle number to device symbol
    try {
	cuda::copy(npart, gpu::ljfluid_square::npart);
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
void ljfluid_gpu_impl_square<dimension>::threads()
{
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
void ljfluid_gpu_impl_square<dimension>::cutoff_radius(float_type value)
{
    try {
	cuda::copy(r_cut, gpu::ljfluid_square::r_cut);
	cuda::copy(rr_cut, gpu::ljfluid_square::rr_cut);
	cuda::copy(en_cut, gpu::ljfluid_square::en_cut);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy potential cutoff symbols");
    }
}

#ifdef USE_POTENTIAL_SMOOTHING
template <int dimension>
void ljfluid_gpu_impl_square<dimension>::potential_smoothing(float_type value)
{
    try {
	cuda::copy(rri_smooth, gpu::ljfluid_square::rri_smooth);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy potential smoothing function scale symbol");
    }
}
#endif /* USE_POTENTIAL_SMOOTHING */

template <int dimension>
void ljfluid_gpu_impl_square<dimension>::box(float_type value)
{
    try {
	cuda::copy(box_, gpu::ljfluid_square::box);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy periodic box length to device symbol");
    }
}

template <int dimension>
void ljfluid_gpu_impl_square<dimension>::timestep(float_type value)
{
    try {
	cuda::copy(timestep_, gpu::ljfluid_square::timestep);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy simulation timestep to device symbol");
    }
}

template <int dimension>
void ljfluid_gpu_impl_square<dimension>::restore(trajectory_visitor visitor)
{
    // read phase space sample
    visitor(m_sample.r, m_sample.v);

    try {
	// copy periodically reduced particle positions from host to GPU
	for (unsigned int i = 0; i < npart; ++i) {
	    h_part.r[i] = make_float(m_sample.r[i]);
	}
	cuda::copy(h_part.r, g_part.r, stream_);
	// replicate to periodically extended particle positions
	cuda::copy(g_part.r, g_part.R, stream_);
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

	// wait for GPU operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to restore system state from phase space sample");
    }
}

template <int dimension>
void ljfluid_gpu_impl_square<dimension>::lattice()
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
void ljfluid_gpu_impl_square<dimension>::temperature(float_type temp)
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

	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to set center of mass velocity to zero");
    }
}

template <int dimension>
void ljfluid_gpu_impl_square<dimension>::mdstep()
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
void ljfluid_gpu_impl_square<dimension>::synchronize()
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
void ljfluid_gpu_impl_square<dimension>::sample()
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
	event_[0].record(stream_);

	// wait for CUDA operations to finish
	event_[0].synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy MD simulation step results from GPU to host");
    }

    std::copy(h_part.r.begin(), h_part.r.end(), m_sample.r.begin());
    std::copy(h_part.R.begin(), h_part.R.end(), m_sample.R.begin());
    std::copy(h_part.v.begin(), h_part.v.end(), m_sample.v.begin());

    // mean potential energy per particle
    m_sample.en_pot = (double) reduce_en.value() / npart;
    // mean virial equation sum per particle
    m_sample.virial = (double) reduce_virial.value() / npart;

    // GPU time for sample memcpy
    m_times["sample_memcpy"] += event_[0] - event_[1];
}

template <int dimension>
void ljfluid_gpu_impl_square<dimension>::attrs(H5::Group const& param) const
{
    // no implementation-dependent parameters
}

template <int dimension>
void ljfluid_gpu_impl_square<dimension>::velocity_verlet(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, stream);
    gpu::ljfluid_square::inteq(g_part.r, g_part.R, g_part.v, g_part.f);
}

template <int dimension>
void ljfluid_gpu_impl_square<dimension>::update_forces(cuda::stream& stream)
{
    cuda::configure(dim_.grid, dim_.block, dim_.threads_per_block() * sizeof(gpu_vector_type), stream);
    gpu::ljfluid_square::mdstep(g_part.r, g_part.v, g_part.f, g_part.en, g_part.virial);
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_GPU_SQUARE_HPP */
