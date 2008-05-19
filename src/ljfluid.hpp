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
#include "exception.hpp"
#include "options.hpp"
#include "rand48.hpp"
#include "vector2d.hpp"
#include "vector3d.hpp"


namespace mdsim
{

/**
 * MD simulation cell placeholders
 */
template <unsigned dimension, typename T>
struct cell_array
{
    /** n-dimensional host floating-point vector type */
    typedef T vector_type;
    /** n-dimensional device floating-point vector type */
    typedef typename cuda::vector_type<dimension, typename T::value_type>::type cuda_vector_type;

    /** n-dimensional particle phase space coordiates */
    cuda::vector<cuda_vector_type> r_gpu, r_gpu2;
    cuda::vector<cuda_vector_type> v_gpu, v_gpu2;
    cuda::host::vector<vector_type> r;
    cuda::host::vector<vector_type> v;
    /** periodically reduced particle coordinates */
    cuda::vector<cuda_vector_type> rp_gpu, rp_gpu2;
    /** particle numbers */
    cuda::vector<int> tag_gpu, tag_gpu2;
    cuda::host::vector<int> tag;
    /** n-dimensional force acting upon particle */
    cuda::vector<cuda_vector_type> force_gpu, force_gpu2;
    cuda::host::vector<vector_type> force;
    /** potential energy and virial equation sum */
    cuda::vector<float2> en_gpu;
    cuda::host::vector<float2> en;

    cell_array() {}

    void resize(size_t n)
    {
	r_gpu.resize(n);
	r_gpu2.resize(n);
	v_gpu.resize(n);
	v_gpu2.resize(n);
	r.resize(n);
	v.resize(n);
	rp_gpu.resize(n);
	rp_gpu2.resize(n);
	tag_gpu.resize(n);
	tag_gpu2.resize(n);
	tag.resize(n);
	force_gpu.resize(n);
	force_gpu2.resize(n);
	force.resize(n);
	en_gpu.resize(n);
	en.resize(n);
    }
};


/**
 * Simulate a Lennard-Jones fluid with naive N-squared algorithm
 */
template <unsigned int NDIM, typename T>
class ljfluid
{
public:
    /** n-dimensional host floating-point vector type */
    typedef T vector_type;
    /** n-dimensional device floating-point vector type */
    typedef typename cuda::vector_type<NDIM, typename T::value_type>::type cuda_vector_type;

public:
    ljfluid(options const& opts);

    uint64_t particles() const;
    float timestep();
    void timestep(float val);
    float density() const;
    void density(float density_);
    float box() const;
    void temperature(float temp);

    void step(float& en_pot, float& virial, T& vel_cm, float& vel2_sum);
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
    /** particle phase space coordinates */
    phase_space_point<cuda::host::vector<T> > part;
    /** CUDA execution dimensions */
    cuda::config dim_;
    /** random number generator */
    mdsim::rand48 rng_;

    /** cell placeholders */
    cell_array<NDIM, T> cell;
    /** CUDA execution dimensions for cell kernels */
    cuda::config cell_dim_;
    /** number of cells per dimension */
    unsigned int ncell_;
    /** total number of cell placeholders */
    unsigned int nplace_;
    /** cell length */
    float r_cell_;

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
};


/**
 * initialize Lennard-Jones fluid with given particle number
 */
template <unsigned int NDIM, typename T>
ljfluid<NDIM, T>::ljfluid(options const& opts) : npart(opts.npart()), part(npart), dim_(opts.dim()), rng_(dim_), steps_(0), gputime_(0.), memtime_(0.)
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

    cuda::copy(npart, gpu::ljfluid::npart);
    cuda::copy(rr_cut, gpu::ljfluid::rr_cut);
    cuda::copy(en_cut, gpu::ljfluid::en_cut);

    // seed random number generator
    rng_.set(opts.rngseed());
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
    cuda::copy(timestep_, gpu::ljfluid::timestep);
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
    cuda::copy(box_, gpu::ljfluid::box);

    // FIXME optimal cell length must consider particle density fluctuations!
    r_cell_ = std::max(r_cut, powf((CELL_SIZE / 4.) / density_, 1.0 / NDIM));
    // optimal number of cells per dimension
    ncell_ = floorf(box_ / r_cell_);
    cuda::copy(ncell_, gpu::ljfluid::ncell);
    // CUDA execution dimensions for cell kernels
    cell_dim_ = cuda::config(dim3(pow(ncell_, NDIM)), dim3(CELL_SIZE));
    // total number of cell placeholders
    nplace_ = pow(ncell_, NDIM) * CELL_SIZE;

    std::cerr << "Placeholders per cell:\t" << CELL_SIZE << "\n";
    std::cerr << "Optimal cell length:\t" << r_cell_ << "\n";
    std::cerr << "Cells per dimension:\t" << ncell_ << "\n";
    std::cerr << "Average cell occupancy:\t" << (float(npart) / nplace_) << std::endl;

    try {
	cell.resize(nplace_);
    }
    catch (cuda::error const& e) {
	throw exception("failed to allocate global device memory for cells");
    }

    // initialize coordinates
    cuda::vector<cuda_vector_type> part(npart);
    gpu::ljfluid::lattice.configure(dim_, stream_);
    gpu::ljfluid::lattice(part.data());

    // assign particles to cells
    gpu::ljfluid::assign_cells.configure(cell_dim_, stream_);
    gpu::ljfluid::assign_cells(part.data(), cell.r_gpu.data(), cell.tag_gpu.data());

    cuda::copy(cell.r_gpu, cell.r, stream_);
    cuda::copy(cell.tag_gpu, cell.tag, stream_);
    cuda::copy(cell.r_gpu, cell.rp_gpu, stream_);
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
void ljfluid<NDIM, T>::temperature(float temp)
{
    // initialize velocities
    cuda::vector<cuda_vector_type> v(npart);
    gpu::ljfluid::boltzmann.configure(dim_, stream_);
    gpu::ljfluid::boltzmann(v.data(), temp, rng_.data());
    cuda::copy(v, part.v, stream_);

    try {
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	throw exception("CUDA ERROR: failed to initialize velocities");
    }

    for (unsigned int i = 0; i < cell.v.size(); ++i) {
	if (IS_REAL_PARTICLE(cell.tag[i])) {
	    cell.v[i] = part.v[cell.tag[i]];
	}
    }
    cuda::copy(cell.v, cell.v_gpu, stream_);

    // initialize forces
    fill(cell.force.begin(), cell.force.end(), 0.);
    cuda::copy(cell.force, cell.force_gpu, stream_);

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

    // integrate equations of motion
    gpu::ljfluid::inteq.configure(cell_dim_, stream_);
    gpu::ljfluid::inteq(cell.r_gpu.data(), cell.rp_gpu.data(), cell.v_gpu.data(), cell.force_gpu.data());

    // update cells
    gpu::ljfluid::update_cells.configure(cell_dim_, stream_);
    gpu::ljfluid::update_cells(cell.r_gpu.data(), cell.rp_gpu.data(), cell.v_gpu.data(), cell.force_gpu.data(), cell.tag_gpu.data(), cell.r_gpu2.data(), cell.rp_gpu2.data(), cell.v_gpu2.data(), cell.force_gpu2.data(), cell.tag_gpu2.data());
    cuda::copy(cell.r_gpu2, cell.r_gpu, stream_);
    cuda::copy(cell.rp_gpu2, cell.rp_gpu, stream_);
    cuda::copy(cell.v_gpu2, cell.v_gpu, stream_);
    cuda::copy(cell.force_gpu2, cell.force_gpu, stream_);
    cuda::copy(cell.tag_gpu2, cell.tag_gpu, stream_);

    // update forces
    gpu::ljfluid::mdstep.configure(cell_dim_, stream_);
    gpu::ljfluid::mdstep(cell.rp_gpu.data(), cell.v_gpu.data(), cell.force_gpu.data(), cell.tag_gpu.data(), cell.en_gpu.data());

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
    cuda::copy(cell.r_gpu, cell.r, stream_);
    cuda::copy(cell.v_gpu, cell.v, stream_);
    cuda::copy(cell.tag_gpu, cell.tag, stream_);
    cuda::copy(cell.en_gpu, cell.en, stream_);
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

    // number of particles found in cells
    unsigned int count = 0;

    for (uint64_t i = 0; i < cell.r.size(); ++i) {
	// check if real particle
	if (IS_REAL_PARTICLE(cell.tag[i])) {
	    // copy phase space coordinates
	    part.r[cell.tag[i]] = cell.r[i];
	    part.v[cell.tag[i]] = cell.v[i];

	    en_pot += (cell.en[i].x - en_pot) / (i + 1);
	    virial += (cell.en[i].y - virial) / (i + 1);
	    vel_cm += (cell.v[i] - vel_cm) / (i + 1);
	    vel2_sum += (cell.v[i] * cell.v[i] - vel2_sum) / (i + 1);

	    count++;
	}
    }

    if (count != npart) throw exception("I've lost my marbles!");
}

/**
 * write particle coordinates and velocities to output stream
 */
template <unsigned int NDIM, typename T>
void ljfluid<NDIM, T>::trajectories(std::ostream& os) const
{
    for (uint64_t i = 0; i < npart; ++i) {
	os << i << "\t" << part.r[i] << "\t" << part.v[i] << "\n";
    }
    os << "\n\n";
}

/**
 * sample trajectory
 */
template <unsigned int NDIM, typename T>
template <typename V>
void ljfluid<NDIM, T>::sample(V& visitor) const
{
    visitor.sample(part);
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
