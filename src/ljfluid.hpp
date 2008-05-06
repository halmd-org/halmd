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
#include "rand48.hpp"
#include "trajectory.hpp"


namespace mdsim
{

/**
 * MD simulation particle
 */
template <typename T>
struct particle
{
    /** n-dimensional particle coordinates */
    cuda::vector<T> pos_gpu;
#ifndef USE_LEAPFROG
    cuda::vector<T> pos_old_gpu;
#endif
    cuda::host::vector<T> pos;
    /** n-dimensional particle velocity */
    cuda::vector<T> vel_gpu;
    cuda::host::vector<T> vel;
    /** n-dimensional force acting upon particle */
    cuda::vector<T> force_gpu;
    cuda::host::vector<T> force;
    /** potential energy */
    cuda::vector<float> en_gpu;
    cuda::host::vector<float> en;
    /** virial equation sum */
    cuda::vector<float> virial_gpu;
    cuda::host::vector<float> virial;

    particle(size_t n) :
	pos_gpu(n),
#ifndef USE_LEAPFROG
	pos_old_gpu(n),
#endif
	pos(n), vel_gpu(n), vel(n), force_gpu(n), force(n),
	en_gpu(n), en(n), virial_gpu(n), virial(n)
    {
    }
};


/**
 * Simulate a Lennard-Jones fluid with naive N-squared algorithm
 */
template <typename T>
class ljfluid
{
public:
    ljfluid(size_t npart, cuda::config const& dim);

    size_t particles() const;
    float timestep();
    void timestep(float val);
    float density() const;
    void density(float density_);
    float box() const;
    void temperature(float temp, rand48& rng);

    void step(float& en_pot, float& virial, T& vel_cm, float& vel2_sum);
    void trajectories(std::ostream& os) const;
    template <unsigned int X, typename Y>
    void trajectories(trajectory<X, Y>& traj) const;

    float gputime() const;
    float memtime() const;

private:
    /** number of particles in periodic box */
    size_t npart;
#ifdef DIM_3D
    /** particles */
    particle<float3> part;
#else
    /** particles */
    particle<float2> part;
#endif
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
};


/**
 * initialize Lennard-Jones fluid with given particle number
 */
template <typename T>
ljfluid<T>::ljfluid(size_t npart, cuda::config const& dim) : npart(npart), part(npart), dim_(dim), steps_(0), gputime_(0.), memtime_(0.)
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
template <typename T>
size_t ljfluid<T>::particles() const
{
    return npart;
}

/**
 * get simulation timestep
 */
template <typename T>
float ljfluid<T>::timestep()
{
    return timestep_;
}

/**
 * set simulation timestep
 */
template <typename T>
void ljfluid<T>::timestep(float timestep_)
{
    this->timestep_ = timestep_;
    gpu::ljfluid::timestep = timestep_;
}

/**
 * get particle density
 */
template <typename T>
float ljfluid<T>::density() const
{
    return density_;
}

/**
 * set particle density
 */
template <typename T>
void ljfluid<T>::density(float density_)
{
    // particle density
    this->density_ = density_;
    // periodic box length
    box_ = pow(npart / density_, 1.0 / T::dim());
    gpu::ljfluid::box = box_;

    // initialize coordinates
    cuda::stream stream;

    gpu::ljfluid::lattice.configure(dim_, stream);
    gpu::ljfluid::lattice(part.pos_gpu.data());

    part.pos.memcpy(part.pos_gpu, stream);
    stream.synchronize();
}

/**
 * get periodic box length
 */
template <typename T>
float ljfluid<T>::box() const
{
    return box_;
}

/**
 * set temperature
 */
template <typename T>
void ljfluid<T>::temperature(float temp, rand48& rng)
{
    cuda::stream stream;

    // initialize velocities
    gpu::ljfluid::boltzmann.configure(dim_, stream);
#ifdef USE_LEAPFROG
    gpu::ljfluid::boltzmann(part.vel_gpu.data(), temp, rng.data());
#else
    gpu::ljfluid::boltzmann(part.vel_gpu.data(), part.pos_gpu.data(), part.pos_old_gpu.data(), temp, rng.data());
#endif
    part.vel.memcpy(part.vel_gpu, stream);

    // initialize forces
#ifdef DIM_3D
    fill(part.force.begin(), part.force.end(), make_float3(0., 0., 0.));
#else
    fill(part.force.begin(), part.force.end(), make_float2(0., 0.));
#endif
    part.force_gpu.memcpy(part.force, stream);

    stream.synchronize();
}

/**
 * MD simulation step
 */
template <typename T>
void ljfluid<T>::step(float& en_pot, float& virial, T& vel_cm, float& vel2_sum)
{
    cuda::stream stream;
    cuda::event start, stop;
    start.record(stream);

    ++steps_;

#ifdef USE_LEAPFROG
    gpu::ljfluid::inteq.configure(dim_, stream);
    gpu::ljfluid::inteq(part.pos_gpu.data(), part.vel_gpu.data(), part.force_gpu.data());
#endif

#ifdef DIM_3D
    // reserve shared device memory for particle coordinates
    gpu::ljfluid::mdstep.configure(dim_, dim_.threads_per_block() * sizeof(float3), stream);
#else
    gpu::ljfluid::mdstep.configure(dim_, dim_.threads_per_block() * sizeof(float2), stream);
#endif
    gpu::ljfluid::mdstep(part.pos_gpu.data(), part.vel_gpu.data(), part.force_gpu.data(), part.en_gpu.data(), part.virial_gpu.data());

#ifndef USE_LEAPFROG
    gpu::ljfluid::inteq.configure(dim_, stream);
    gpu::ljfluid::inteq(part.pos_gpu.data(), part.pos_old_gpu.data(), part.vel_gpu.data(), part.force_gpu.data());
#endif

    stop.record(stream);
    stop.synchronize();

    // adjust average GPU time in milliseconds per simulation time
    gputime_ += ((stop - start) - gputime_) / steps_;

    start.record(stream);

    part.pos.memcpy(part.pos_gpu, stream);
    part.vel.memcpy(part.vel_gpu, stream);
    part.force.memcpy(part.force_gpu, stream);
    part.en.memcpy(part.en_gpu, stream);
    part.virial.memcpy(part.virial_gpu, stream);

    stop.record(stream);
    stop.synchronize();

    // adjust average device memory transfer time in milliseconds per simulation time
    memtime_ += ((stop - start) - memtime_) / steps_;

    // compute averages
    en_pot = 0.;
    virial = 0.;
    vel_cm = 0.;
    vel2_sum = 0.;

    for (size_t i = 0; i < npart; ++i) {
	en_pot += (part.en[i] - en_pot) / (i + 1);
	virial += (part.virial[i] - virial) / (i + 1);
	vel_cm += (T(part.vel[i]) - vel_cm) / (i + 1);
	vel2_sum += (T(part.vel[i]) * T(part.vel[i]) - vel2_sum) / (i + 1);
    }
}

/**
 * write particle coordinates and velocities to output stream
 */
template <typename T>
void ljfluid<T>::trajectories(std::ostream& os) const
{
    for (size_t i = 0; i < npart; ++i) {
	os << T(part.pos[i]) << "\t" << T(part.vel[i]) << "\n";
    }
    os << "\n\n";
}

/**
 * write particle coordinates and velocities to binary HDF5 file
 */
template <typename T>
template <unsigned int X, typename Y>
void ljfluid<T>::trajectories(trajectory<X, Y>& traj) const
{
    traj.write(part.pos, part.vel, part.force);
}

/**
 * get total GPU time in seconds
 */
template <typename T>
float ljfluid<T>::gputime() const
{
    return (gputime_ * 1.e-3) * steps_;
}

/**
 * get total device memory transfer time in seconds
 */
template <typename T>
float ljfluid<T>::memtime() const
{
    return (memtime_ * 1.e-3) * steps_;
}

} // namespace mdsim

#endif /* ! MDSIM_LJFLUID_HPP */
