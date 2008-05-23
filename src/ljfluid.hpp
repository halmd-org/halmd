/* Lennard-Jones fluid
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
#include <cuda_wrapper.hpp>
#include <math.h>
#include <stdint.h>
#include "gpu/ljfluid_glue.hpp"
#include "exception.hpp"
#include "options.hpp"
#include "statistics.hpp"
#include "rand48.hpp"
#include "vector2d.hpp"
#include "vector3d.hpp"


namespace mdsim
{

/**
 * MD simulation particle
 */
template <unsigned dimension, typename T>
struct particle
{
    /** device floating-point vector type */
    typedef typename cuda::types::vector<dimension, typename T::value_type>::type floatn;

    /** n-dimensional particle phase space coordiates */
    phase_space_point<cuda::vector<floatn> > psc_gpu;
    phase_space_point<cuda::host::vector<T> > psc;
    /** periodically reduced particle coordinates */
    cuda::vector<floatn> rp_gpu;
    /** n-dimensional force acting upon particle */
    cuda::vector<floatn> force_gpu;
    cuda::host::vector<T> force;
    /** potential energy and virial equation sum */
    cuda::vector<float2> en_gpu;
    cuda::host::vector<vector2d<float> > en;

    particle(uint64_t n) : psc_gpu(n), psc(n), rp_gpu(n), force_gpu(n), force(n), en_gpu(n), en(n) { }
};


/**
 * Simulate a Lennard-Jones fluid with naive N-squared algorithm
 */
template <unsigned dimension, typename T>
class ljfluid
{
public:
    /** device floating-point vector type */
    typedef typename cuda::types::vector<dimension, typename T::value_type>::type floatn;

public:
    ljfluid(options const& opts);

    uint64_t particles() const;
    float timestep();
    void timestep(float val);
    float density() const;
    void density(float density_);
    float box() const;
    void temperature(float temp);

    void mdstep();
    void trajectories(std::ostream& os) const;
    template <typename V>
    void sample(V& visitor) const;

    float gputime() const;
    float memtime() const;

private:
    void step();

private:
    /** number of particles in periodic box */
    uint64_t npart;
    /** particles */
    particle<dimension, T> part;
    /** random number generator */
    mdsim::rand48 rng_;
    /** CUDA execution dimensions */
    cuda::config dim_;

    /** particle density */
    float density_;
    /** periodic box length */
    float box_;
    /** MD simulation timestep */
    float timestep_;
    /** cutoff distance for shifted Lennard-Jones potential */
    float r_cut;

    /** total number of simulation steps */
    uint64_t steps_;
    /** average GPU time in milliseconds per simulation step */
    float gputime_;
    /** average device memory transfer time in milliseconds per simulation step */
    float memtime_;

    /** CUDA asynchronous execution */
    cuda::stream stream_;
    cuda::event event_[4];

    /** potential energy per particle */
    float en_pot_;
    /** virial theorem force sum */
    float virial_;
};


/**
 * initialize Lennard-Jones fluid with given particle number
 */
template <unsigned dimension, typename T>
ljfluid<dimension, T>::ljfluid(options const& opts) : npart(opts.npart()), part(opts.npart()), rng_(opts.dim()), dim_(opts.dim()), steps_(0), gputime_(0.), memtime_(0.)
{
    // fixed cutoff distance for shifted Lennard-Jones potential
    r_cut = 2.5;
    // squared cutoff distance
    float rr_cut = r_cut * r_cut;

    // potential energy at cutoff distance
    float rri_cut = 1. / rr_cut;
    float r6i_cut = rri_cut * rri_cut * rri_cut;
    float en_cut = 4. * r6i_cut * (r6i_cut - 1.);

    cuda::copy(npart, gpu::ljfluid::npart);
    cuda::copy(rr_cut, gpu::ljfluid::rr_cut);
    cuda::copy(en_cut, gpu::ljfluid::en_cut);

    // seed random number generator
    rng_.set(opts.rngseed());

    // reserve device memory for placeholder particles
    part.psc_gpu.r.reserve(dim_.threads());
    part.psc_gpu.v.reserve(dim_.threads());
    part.rp_gpu.reserve(dim_.threads());
    part.force_gpu.reserve(dim_.threads());
    part.en_gpu.reserve(dim_.threads());
}

/**
 * get number of particles in periodic box
 */
template <unsigned dimension, typename T>
uint64_t ljfluid<dimension, T>::particles() const
{
    return npart;
}

/**
 * get simulation timestep
 */
template <unsigned dimension, typename T>
float ljfluid<dimension, T>::timestep()
{
    return timestep_;
}

/**
 * set simulation timestep
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::timestep(float timestep_)
{
    this->timestep_ = timestep_;
    cuda::copy(timestep_, gpu::ljfluid::timestep);
}

/**
 * get particle density
 */
template <unsigned dimension, typename T>
float ljfluid<dimension, T>::density() const
{
    return density_;
}

/**
 * set particle density
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::density(float density_)
{
    // particle density
    this->density_ = density_;
    // periodic box length
    box_ = pow(npart / density_, 1.0 / dimension);
    cuda::copy(box_, gpu::ljfluid::box);

    // initialize coordinates
    gpu::ljfluid::lattice.configure(dim_, stream_);
    gpu::ljfluid::lattice(part.psc_gpu.r.data());

    cuda::copy(part.psc_gpu.r, part.psc.r, stream_);
    cuda::copy(part.psc_gpu.r, part.rp_gpu, stream_);
    stream_.synchronize();
}

/**
 * get periodic box length
 */
template <unsigned dimension, typename T>
float ljfluid<dimension, T>::box() const
{
    return box_;
}

/**
 * set temperature
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::temperature(float temp)
{
    // initialize velocities
    gpu::ljfluid::boltzmann.configure(dim_, stream_);
    gpu::ljfluid::boltzmann(part.psc_gpu.v.data(), temp, rng_.data());
    cuda::copy(part.psc_gpu.v, part.psc.v, stream_);

    try {
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("CUDA ERROR: failed to initialize velocities");
    }

    // set center of mass velocity to zero
    T v_cm = mean(part.psc.v.begin(), part.psc.v.end());
    for (typename cuda::host::vector<T>::iterator v = part.psc.v.begin(); v != part.psc.v.end(); ++v) {
	*v -= v_cm;
    }
    cuda::copy(part.psc.v, part.psc_gpu.v, stream_);

    // initialize forces
    fill(part.force.begin(), part.force.end(), 0.);
    cuda::copy(part.force, part.force_gpu, stream_);

    try {
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("CUDA ERROR: failed to initialize forces");
    }
}

/**
 * MD simulation step on GPU
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::step()
{
    event_[0].record(stream_);
    gpu::ljfluid::inteq.configure(dim_, stream_);
    gpu::ljfluid::inteq(part.psc_gpu.r.data(), part.rp_gpu.data(), part.psc_gpu.v.data(), part.force_gpu.data());
    // reserve shared device memory for particle coordinates
    gpu::ljfluid::mdstep.configure(dim_, dim_.threads_per_block() * sizeof(T), stream_);
    gpu::ljfluid::mdstep(part.rp_gpu.data(), part.psc_gpu.v.data(), part.force_gpu.data(), part.en_gpu.data());
    event_[1].record(stream_);
}

/**
 * MD simulation step
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::mdstep()
{
    if (steps_ == 0) {
	// first MD simulation step on GPU
	step();
    }

    event_[2].record(stream_);
    cuda::copy(part.psc_gpu.r, part.psc.r, stream_);
    cuda::copy(part.psc_gpu.v, part.psc.v, stream_);
    cuda::copy(part.en_gpu, part.en, stream_);
    event_[3].record(stream_);
    event_[3].synchronize();

    steps_++;

    // adjust average GPU time in milliseconds per simulation time
    gputime_ += ((event_[1] - event_[0]) - gputime_) / steps_;
    // adjust average device memory transfer time in milliseconds per simulation time
    memtime_ += ((event_[3] - event_[2]) - memtime_) / steps_;

    // next MD simulation step on GPU
    step();

    // compute average potential energy and virial theorem sum per particle
    vector2d<typename T::value_type> en = mean(part.en.begin(), part.en.end());
    // average potential energy per particle
    en_pot_ = en.x;
    // average virial theorem sum per particle
    virial_ = en.y;
}

/**
 * write particle coordinates and velocities to output stream
 */
template <unsigned dimension, typename T>
void ljfluid<dimension, T>::trajectories(std::ostream& os) const
{
    for (uint64_t i = 0; i < npart; ++i) {
	os << i << "\t" << part.psc.r[i] << "\t" << part.psc.v[i] << "\n";
    }
    os << "\n\n";
}

/**
 * sample trajectory
 */
template <unsigned dimension, typename T>
template <typename V>
void ljfluid<dimension, T>::sample(V& visitor) const
{
    visitor.sample(part.psc, en_pot_, virial_);
}

/**
 * get total GPU time in seconds
 */
template <unsigned dimension, typename T>
float ljfluid<dimension, T>::gputime() const
{
    return (gputime_ * 1.e-3) * steps_;
}

/**
 * get total device memory transfer time in seconds
 */
template <unsigned dimension, typename T>
float ljfluid<dimension, T>::memtime() const
{
    return (memtime_ * 1.e-3) * steps_;
}

} // namespace mdsim

#endif /* ! MDSIM_LJFLUID_HPP */
