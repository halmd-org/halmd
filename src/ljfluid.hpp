/* Lennard-Jones fluid simulation with naive N-squared algorithm
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

#ifndef MDSIM_LJFLUID_HPP
#define MDSIM_LJFLUID_HPP

#include <algorithm>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <cmath>
#include <cuda_wrapper.hpp>
#include <stdint.h>
#include "gpu/ljfluid_glue.hpp"
#include "exception.hpp"
#include "log.hpp"
#include "options.hpp"
#include "statistics.hpp"
#include "rand48.hpp"
#include "vector2d.hpp"
#include "vector3d.hpp"


#define foreach BOOST_FOREACH

namespace mdsim
{

/**
 * Lennard-Jones fluid simulation with naive N-squared algorithm
 */
template <unsigned dimension, typename T>
class ljfluid
{
public:
    /** device floating-point vector type */
    typedef typename cuda::types::vector<dimension, typename T::value_type>::type floatn;

public:
    /** initialize fixed simulation parameters */
    ljfluid();
    /** set number of particles in system */
    void particles(unsigned int value);
    /** set system state from phase space sample */
    void particles(phase_space_point<std::vector<T> > const& state);
    /** set number of CUDA execution threads */
    void threads(unsigned int value);
    /** initialize random number generator with seed */
    void rng(unsigned int seed);
    /** initialize random number generator from state */
    void rng(mdsim::rand48::state_type const& state);

    /** set particle density */
    void density(float value);
    /** set periodic box length */
    void box(float value);
    /** arrange particles on a face-centered cubic (fcc) lattice */
    void lattice();
    /** set system temperature according to Maxwell-Boltzmann distribution */
    void temperature(float temp);
    /** set simulation timestep */
    void timestep(float value);

    /** get number of particles */
    unsigned int const& particles() const;
    /** get particle density */
    float const& density() const;
    /** get periodic box length */
    float const& box() const;
    /** get simulation timestep */
    float const& timestep() const;

    /** stream MD simulation step on GPU */
    void mdstep();
    /** copy MD simulation step results from GPU to host */
    void sample();
    /** sample trajectory */
    template <typename visitor>
    void sample(visitor& v) const;

private:
    /** number of particles in system */
    unsigned int npart;
    /** particle density */
    float density_;
    /** periodic box length */
    float box_;
    /** simulation timestep */
    float timestep_;
    /** cutoff distance for shifted Lennard-Jones potential */
    float r_cut;

    /** system state in phase space */
    phase_space_point<cuda::host::vector<T> > h_state;
    phase_space_point<cuda::vector<floatn> > g_state;
    /** forces */
    cuda::vector<floatn> g_force;
    /** potential energy and virial equation sum */
    cuda::host::vector<vector2d<float> > h_en;
    cuda::vector<float2> g_en;

    /** random number generator */
    mdsim::rand48 rng_;
    /** CUDA execution dimensions */
    cuda::config dim_;
    /** CUDA asynchronous execution */
    cuda::stream stream_;
};


/**
 * initialize fixed simulation parameters
 */
template <unsigned dimension, typename T>
ljfluid<dimension, T>::ljfluid()
{
    // fixed cutoff distance for shifted Lennard-Jones potential
    r_cut = 2.5;

    // squared cutoff distance
    float rr_cut = r_cut * r_cut;
    // potential energy at cutoff distance
    float rri_cut = 1. / rr_cut;
    float r6i_cut = rri_cut * rri_cut * rri_cut;
    float en_cut = 4. * r6i_cut * (r6i_cut - 1.);

    try {
	cuda::copy(rr_cut, gpu::ljfluid::rr_cut);
	cuda::copy(en_cut, gpu::ljfluid::en_cut);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy cutoff parameters to device symbols");
    }
}

/**
 * set number of particles in system
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::particles(unsigned int value)
{
    // validate particle number
    if (value < 1) {
	throw exception("invalid number of particles");
    }
    // set particle number
    npart = value;
    // copy particle number to device symbol
    try {
	cuda::copy(npart, gpu::ljfluid::npart);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy particle number to device symbol");
    }

    // allocate global device memory for system state
    try {
	g_state.r.resize(npart);
	g_state.R.resize(npart);
	g_state.v.resize(npart);
	g_force.resize(npart);
	g_en.resize(npart);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for system state");
    }

    // allocate page-locked host memory for system state
    try {
	h_state.r.resize(npart);
	h_state.R.resize(npart);
	h_state.v.resize(npart);
	// particle forces reside only in GPU memory
	h_en.resize(npart);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate page-locked host memory for system state");
    }

    // set particle forces to zero for first integration of differential equations of motion
    try {
	fill(h_state.v.begin(), h_state.v.end(), 0.);
	cuda::copy(h_state.v, g_force, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to initialize forces to zero");
    }
}

/**
 * set system state from phase space sample
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::particles(phase_space_point<std::vector<T> > const& state)
{
    // set number of particles in system
    particles(state.r.size());

    // set system state from phase space sample
    std::copy(state.r.begin(), state.r.end(), h_state.r.begin());
    try {
	// copy periodically reduced particles positions from host to GPU
	cuda::copy(h_state.r, g_state.r, stream_);
	// replicate to periodically extended particle positions
	cuda::copy(g_state.r, g_state.R, stream_);
	cuda::copy(h_state.r, h_state.R, stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to set system state from phase space sample");
    }
}

/**
 * set number of CUDA execution threads
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::threads(unsigned int value)
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

    // allocate additional global device memory for system state
    try {
	g_state.r.reserve(npart);
	g_state.R.reserve(npart);
	g_state.v.reserve(npart);
	g_force.reserve(npart);
	g_en.reserve(npart);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate additional global device memory for system state");
    }
}

/**
 * initialize random number generator with seed
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::rng(unsigned int seed)
{
    // initialize random number generator
    try {
	rng_.resize(dim_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for random number generator");
    }

    // seed random number generator
    try {
	rng_.set(seed);
    }
    catch (cuda::error const& e) {
	throw exception("failed to seed random number generator");
    }
}

/**
 * initialize random number generator from state
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::rng(mdsim::rand48::state_type const& state)
{
    // initialize random number generator
    try {
	rng_.resize(dim_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for random number generator");
    }

    // restore random number generator state
    try {
	rng_.restore(state);
    }
    catch (cuda::error const& e) {
	throw exception("failed to restore random number generator state");
    }
}

/**
 * set particle density
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::density(float value)
{
    // set particle density
    density_ = value;

    // compute periodic box length
    box_ = powf(npart / density_, 1. / dimension);
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
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::box(float value)
{
    // set periodic box length
    box_ = value;
    // copy periodic box length to device symbol
    try {
	cuda::copy(box_, gpu::ljfluid::box);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy periodic box length to device symbol");
    }

    // compute particle density
    density_ = npart / powf(box_, dimension);
}

/**
 * arrange particles on a face-centered cubic (fcc) lattice
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::lattice()
{
    try {
	// compute particle lattice positions on GPU
	gpu::ljfluid::lattice.configure(dim_, stream_);
	gpu::ljfluid::lattice(g_state.r.data());
	// copy particle positions from GPU to host
	cuda::copy(g_state.r, h_state.r, stream_);
	// copy particle positions to periodically extended positions
	cuda::copy(g_state.r, g_state.R, stream_);
	cuda::copy(h_state.r, h_state.R, stream_);

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute particle lattice positions on GPU");
    }
}

/**
 * set system temperature according to Maxwell-Boltzmann distribution
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::temperature(float temp)
{
    try {
	// set velocities using Maxwell-Boltzmann distribution at temperature
	gpu::ljfluid::boltzmann.configure(dim_, stream_);
	gpu::ljfluid::boltzmann(g_state.v.data(), temp, rng_.data());
	// copy particle velocities from GPU to host
	cuda::copy(g_state.v, h_state.v, stream_);

	// wait for CUDA operations to finish
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to compute Maxwell-Boltzmann distributed velocities on GPU");
    }

    // compute center of mass velocity
    T v_cm = mean(h_state.v.begin(), h_state.v.end());
    // set center of mass velocity to zero
    foreach (T& v, h_state.v) {
	v -= v_cm;
    }
    // copy particle velocities from host to GPU
    try {
	cuda::copy(h_state.v, g_state.v, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("failed to set center of mass velocity to zero");
    }
}

/**
 * set simulation timestep
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::timestep(float value)
{
    // set simulation timestep
    timestep_ = value;
    // copy simulation timestep to device symbol
    try {
	cuda::copy(timestep_, gpu::ljfluid::timestep);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy simulation timestep to device symbol");
    }
}

/**
 * get number of particles
 */
template <unsigned dimension, typename T>
unsigned int const& ljfluid<dimension, T>::particles() const
{
    return npart;
}

/**
 * get particle density
 */
template <unsigned dimension, typename T>
float const& ljfluid<dimension, T>::density() const
{
    return density_;
}

/**
 * get periodic box length
 */
template <unsigned dimension, typename T>
float const& ljfluid<dimension, T>::box() const
{
    return box_;
}

/**
 * get simulation timestep
 */
template <unsigned dimension, typename T>
float const& ljfluid<dimension, T>::timestep() const
{
    return timestep_;
}

/**
 * stream MD simulation step on GPU
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::mdstep()
{
    // first leapfrog step of integration of differential equations of motion
    try {
	gpu::ljfluid::inteq.configure(dim_, stream_);
	gpu::ljfluid::inteq(g_state.r.data(), g_state.R.data(), g_state.v.data(), g_force.data());
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream first leapfrog step on GPU");
    }

    // Lennard-Jones force calculation
    try {
	gpu::ljfluid::mdstep.configure(dim_, dim_.threads_per_block() * sizeof(T), stream_);
	gpu::ljfluid::mdstep(g_state.r.data(), g_state.v.data(), g_force.data(), g_en.data());
    }
    catch (cuda::error const& e) {
	throw exception("failed to stream force calculation on GPU");
    }
}

/**
 * copy MD simulation step results from GPU to host
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::sample()
{
    // wait for MD simulation step on GPU to finish
    try {
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("MD simulation step on GPU failed");
    }

    // copy MD simulation step results from GPU to host
    try {
	// copy periodically reduce particles positions
	cuda::copy(g_state.r, h_state.r, stream_);
	// copy periodically extended particles positions
	cuda::copy(g_state.R, h_state.R, stream_);
	// copy particle velocities
	cuda::copy(g_state.v, h_state.v, stream_);
	// copy potential energies and virial equation sums
	cuda::copy(g_en, h_en, stream_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to copy MD simulation step results from GPU to host");
    }
}

/**
 * sample trajectory
 */
template <unsigned dimension, typename T>
template <typename visitor>
void ljfluid<dimension, T>::sample(visitor& v) const
{
    // compute average potential energy and virial theorem sum per particle
    vector2d<float> en = mean(h_en.begin(), h_en.end());

    // sample trajectory, potential energy and virial equation sum
    v.sample(h_state, en.x, en.y);
}

} // namespace mdsim

#undef foreach

#endif /* ! MDSIM_LJFLUID_HPP */
