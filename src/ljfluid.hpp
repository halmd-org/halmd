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

#include <math.h>
#include <stdint.h>
#include <algorithm>
#include <cuda_wrapper.hpp>
#include "gpu/ljfluid_glue.hpp"
#include "autocorrelation.hpp"
#include "exception.hpp"
#include "rand48.hpp"
#include "trajectory.hpp"
#include "vector2d.hpp"
#include "vector3d.hpp"


namespace mdsim
{

/**
 * MD simulation particle
 */
template <typename T>
struct particle
{
    /** n-dimensional particle phase space coordiates */
    phase_space_point<cuda::vector<T> > psc_gpu;
    phase_space_point<cuda::host::vector<T> > psc;
    /** periodically reduced particle coordinates */
    cuda::vector<T> rp_gpu;
    /** n-dimensional force acting upon particle */
    cuda::vector<T> force_gpu;
    cuda::host::vector<T> force;
    /** potential energy */
    cuda::vector<float> en_gpu;
    cuda::host::vector<float> en;
    /** virial equation sum */
    cuda::vector<float> virial_gpu;
    cuda::host::vector<float> virial;

    particle(uint64_t n) : psc_gpu(n), psc(n), rp_gpu(n), force_gpu(n), force(n), en_gpu(n), en(n), virial_gpu(n), virial(n) { }
};


/**
 * Simulate a Lennard-Jones fluid with naive N-squared algorithm
 */
template <unsigned int NDIM, typename T>
class ljfluid
{
public:
    ljfluid(uint64_t npart, cuda::config const& dim);

    uint64_t particles() const;
    float timestep();
    void timestep(float val);
    float density() const;
    void density(float density_);
    float box() const;
    void temperature(float temp, rand48& rng);

    void step(float& en_pot, float& virial, T& vel_cm, float& vel2_sum);
    void trajectories(std::ostream& os) const;
    template <typename Y>
    void trajectories(trajectory<NDIM, Y>& traj) const;
    void sample(autocorrelation<NDIM, T>& tcf) const;

    float gputime() const;
    float memtime() const;

private:
    void step();

private:
    /** number of particles in periodic box */
    uint64_t npart;
    /** particles */
    particle<T> part;
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

    cuda::stream stream_;
    cuda::event event_[4];
};


/**
 * initialize Lennard-Jones fluid with given particle number
 */
template <unsigned int NDIM, typename T>
ljfluid<NDIM, T>::ljfluid(uint64_t npart, cuda::config const& dim) : npart(npart), part(npart), dim_(dim), steps_(0), gputime_(0.), memtime_(0.)
{
    // FIXME do without this requirement
    assert(npart == dim_.threads());

    // fixed cutoff distance for shifted Lennard-Jones potential
    r_cut = 2.5;
    // squared cutoff distance
    float rr_cut = r_cut * r_cut;

    // potential energy at cutoff distance
    float rri_cut = 1. / rr_cut;
    float r6i_cut = rri_cut * rri_cut * rri_cut;
    float en_cut = 4. * r6i_cut * (r6i_cut - 1.);

    gpu::ljfluid::npart = npart;
    gpu::ljfluid::rr_cut = rr_cut;
    gpu::ljfluid::en_cut = en_cut;
}

/**
 * get number of particles in periodic box
 */
template <unsigned int NDIM, typename T>
uint64_t ljfluid<NDIM, T>::particles() const
{
    return npart;
}

/**
 * get simulation timestep
 */
template <unsigned int NDIM, typename T>
float ljfluid<NDIM, T>::timestep()
{
    return timestep_;
}

/**
 * set simulation timestep
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::timestep(float timestep_)
{
    this->timestep_ = timestep_;
    gpu::ljfluid::timestep = timestep_;
}

/**
 * get particle density
 */
template <unsigned int NDIM, typename T>
float ljfluid<NDIM, T>::density() const
{
    return density_;
}

/**
 * set particle density
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::density(float density_)
{
    // particle density
    this->density_ = density_;
    // periodic box length
    box_ = pow(npart / density_, 1.0 / NDIM);
    gpu::ljfluid::box = box_;

    // initialize coordinates
    gpu::ljfluid::lattice.configure(dim_, stream_);
    gpu::ljfluid::lattice(cuda_cast(part.psc_gpu.r));

    part.psc.r.memcpy(part.psc_gpu.r, stream_);
    part.rp_gpu.memcpy(part.psc_gpu.r, stream_);
    stream_.synchronize();
}

/**
 * get periodic box length
 */
template <unsigned int NDIM, typename T>
float ljfluid<NDIM, T>::box() const
{
    return box_;
}

/**
 * set temperature
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::temperature(float temp, rand48& rng)
{
    // initialize velocities
    gpu::ljfluid::boltzmann.configure(dim_, stream_);
    gpu::ljfluid::boltzmann(cuda_cast(part.psc_gpu.v), temp, cuda_cast(rng));
    part.psc.v.memcpy(part.psc_gpu.v, stream_);

    try {
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("CUDA ERROR: failed to initialize velocities");
    }

    // initialize forces
    fill(part.force.begin(), part.force.end(), 0.);
    part.force_gpu.memcpy(part.force, stream_);

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
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::step()
{
    event_[0].record(stream_);
    gpu::ljfluid::inteq.configure(dim_, stream_);
    gpu::ljfluid::inteq(cuda_cast(part.psc_gpu.r), cuda_cast(part.rp_gpu), cuda_cast(part.psc_gpu.v), cuda_cast(part.force_gpu));
    // reserve shared device memory for particle coordinates
    gpu::ljfluid::mdstep.configure(dim_, dim_.threads_per_block() * sizeof(T), stream_);
    gpu::ljfluid::mdstep(cuda_cast(part.rp_gpu), cuda_cast(part.psc_gpu.v), cuda_cast(part.force_gpu), cuda_cast(part.en_gpu), cuda_cast(part.virial_gpu));
    event_[1].record(stream_);
}

/**
 * MD simulation step
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::step(float& en_pot, float& virial, T& vel_cm, float& vel2_sum)
{
    if (steps_ == 0) {
	// first MD simulation step on GPU
	step();
    }

    event_[2].record(stream_);
    part.psc.r.memcpy(part.psc_gpu.r, stream_);
    part.psc.v.memcpy(part.psc_gpu.v, stream_);
    part.en.memcpy(part.en_gpu, stream_);
    part.virial.memcpy(part.virial_gpu, stream_);
    event_[3].record(stream_);
    event_[3].synchronize();

    steps_++;

    // adjust average GPU time in milliseconds per simulation time
    gputime_ += ((event_[1] - event_[0]) - gputime_) / steps_;
    // adjust average device memory transfer time in milliseconds per simulation time
    memtime_ += ((event_[3] - event_[2]) - memtime_) / steps_;

    // next MD simulation step on GPU
    step();

    // compute averages
    en_pot = 0.;
    virial = 0.;
    vel_cm = 0.;
    vel2_sum = 0.;

    for (uint64_t i = 0; i < npart; ++i) {
	en_pot += (part.en[i] - en_pot) / (i + 1);
	virial += (part.virial[i] - virial) / (i + 1);
	vel_cm += (part.psc.v[i] - vel_cm) / (i + 1);
	vel2_sum += (part.psc.v[i] * part.psc.v[i] - vel2_sum) / (i + 1);
    }
}

/**
 * write particle coordinates and velocities to output stream
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::trajectories(std::ostream& os) const
{
    for (uint64_t i = 0; i < npart; ++i) {
	os << part.psc.r[i] << "\t" << part.psc.v[i] << "\n";
    }
    os << "\n\n";
}

/**
 * write particle coordinates and velocities to binary HDF5 file
 */
template <unsigned int NDIM, typename T>
template <typename Y>
void ljfluid<NDIM, T>::trajectories(trajectory<NDIM, Y>& traj) const
{
    traj.write(part.psc);
}

/**
 * sample trajectory for autocorrelation
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::sample(autocorrelation<NDIM, T>& tcf) const
{
    tcf.sample(part.psc);
}

/**
 * get total GPU time in seconds
 */
template <unsigned int NDIM, typename T>
float ljfluid<NDIM, T>::gputime() const
{
    return (gputime_ * 1.e-3) * steps_;
}

/**
 * get total device memory transfer time in seconds
 */
template <unsigned int NDIM, typename T>
float ljfluid<NDIM, T>::memtime() const
{
    return (memtime_ * 1.e-3) * steps_;
}

} // namespace mdsim

#endif /* ! MDSIM_LJFLUID_HPP */
