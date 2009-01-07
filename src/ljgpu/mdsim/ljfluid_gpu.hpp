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

#ifndef LJGPU_MDSIM_LJFLUID_GPU_HPP
#define LJGPU_MDSIM_LJFLUID_GPU_HPP

#include <ljgpu/mdsim/ljfluid_gpu_square.hpp>
#include <ljgpu/mdsim/ljfluid_gpu_cell.hpp>
#include <ljgpu/mdsim/ljfluid_gpu_nbr.hpp>

namespace ljgpu
{

template <template <int> class ljfluid_gpu_impl, int dimension>
class ljfluid_gpu : public ljfluid_gpu_impl<dimension>
{
public:
    typedef typename ljfluid_gpu_impl<dimension>::float_type float_type;
    typedef typename ljfluid_gpu_impl<dimension>::vector_type vector_type;
    typedef typename ljfluid_gpu_impl<dimension>::gpu_vector_type gpu_vector_type;

public:
    /** set potential cutoff radius */
    void cutoff_radius(float_type value);
#ifdef USE_POTENTIAL_SMOOTHING
    /** set potential smoothing function scale parameter */
    void potential_smoothing(float_type value);
#endif
    /** set number of particles in system */
    void particles(unsigned int value);
    /** set number of CUDA execution threads */
    void threads(unsigned int value);
    /* set particle density */
    void density(float_type value);
    /* set periodic box length */
    void box(float_type value);
    /** set simulation timestep */
    void timestep(float_type value);

    /** seed random number generator */
    void rng(unsigned int seed);
    /** restore random number generator from state */
    void rng(rand48::state_type const& state);

    /** get potential cutoff radius */
    float_type cutoff_radius() const { return r_cut; }
#ifdef USE_POTENTIAL_SMOOTHING
    /** get potential smoothing function scale parameter */
    float_type potential_smoothing() const { return r_smooth; }
#endif
    /** get number of particles */
    unsigned int particles() const { return npart; }
    /** get number of CUDA execution blocks */
    unsigned int blocks() const { return dim_.blocks_per_grid(); }
    /** get number of CUDA execution threads */
    unsigned int threads() const { return dim_.threads_per_block(); }
    /** get particle density */
    float_type density() const { return density_; }
    /** get periodic box length */
    float_type box() const { return box_; }
    /** get simulation timestep */
    float_type timestep() const { return timestep_; }
    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;

protected:
    using ljfluid_gpu_impl<dimension>::r_cut;
#ifdef USE_POTENTIAL_SMOOTHING
    using ljfluid_gpu_impl<dimension>::r_smooth;
    using ljfluid_gpu_impl<dimension>::rri_smooth;
#endif
    using ljfluid_gpu_impl<dimension>::rr_cut;
    using ljfluid_gpu_impl<dimension>::en_cut;
    using ljfluid_gpu_impl<dimension>::npart;
    using ljfluid_gpu_impl<dimension>::density_;
    using ljfluid_gpu_impl<dimension>::box_;
    using ljfluid_gpu_impl<dimension>::timestep_;
    using ljfluid_gpu_impl<dimension>::rng_;
    using ljfluid_gpu_impl<dimension>::m_sample;
    using ljfluid_gpu_impl<dimension>::dim_;
    using ljfluid_gpu_impl<dimension>::stream_;
};

template <template <int> class ljfluid_gpu_impl, int dimension>
void ljfluid_gpu<ljfluid_gpu_impl, dimension>::cutoff_radius(float_type value)
{
    r_cut = value;
    LOG("potential cutoff radius: " << r_cut);

    // squared cutoff radius
    rr_cut = std::pow(r_cut, 2);
    // potential energy at cutoff radius
    float_type rri_cut = 1 / rr_cut;
    float_type r6i_cut = rri_cut * rri_cut * rri_cut;
    en_cut = 4 * r6i_cut * (r6i_cut - 1);

    LOG("potential cutoff energy: " << en_cut);

    ljfluid_gpu_impl<dimension>::cutoff_radius(r_cut);
}

#ifdef USE_POTENTIAL_SMOOTHING
template <template <int> class ljfluid_gpu_impl, int dimension>
void ljfluid_gpu<ljfluid_gpu_impl, dimension>::potential_smoothing(float_type value)
{
    r_smooth = value;
    LOG("potential smoothing function scale parameter: " << r_smooth);

    // squared inverse potential smoothing function scale parameter
    rri_smooth = std::pow(r_smooth, -2);

    ljfluid_gpu_impl<dimension>::potential_smoothing(r_cut);
}
#endif /* USE_POTENTIAL_SMOOTHING */

template <template <int> class ljfluid_gpu_impl, int dimension>
void ljfluid_gpu<ljfluid_gpu_impl, dimension>::particles(unsigned int value)
{
    // validate particle number
    if (value < 1) {
	throw exception("invalid number of particles");
    }
    // set particle number
    npart = value;
    LOG("number of particles: " << npart);

    // allocate swappable host memory for trajectory sample
    try {
	m_sample.r.resize(npart);
	m_sample.R.resize(npart);
	m_sample.v.resize(npart);
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate swappable host memory for trajectory sample");
    }

    // implementation-dependent memory allocation
    ljfluid_gpu_impl<dimension>::particles(npart);
}

template <template <int> class ljfluid_gpu_impl, int dimension>
void ljfluid_gpu<ljfluid_gpu_impl, dimension>::threads(unsigned int value)
{
    // query CUDA device properties
    cuda::device::properties prop;
    try {
	prop = cuda::device::properties(cuda::device::get());
    }
    catch (cuda::error const&) {
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
    LOG("number of CUDA execution blocks: " << dim_.blocks_per_grid());
    LOG("number of CUDA execution threads: " << dim_.threads_per_block());

    // change random number generator dimensions
    try {
	rng_.resize(dim_);
    }
    catch (cuda::error const&) {
	throw exception("failed to change random number generator dimensions");
    }

    // implementation-dependent thread allocation
    ljfluid_gpu_impl<dimension>::threads();
}

template <template <int> class ljfluid_gpu_impl, int dimension>
void ljfluid_gpu<ljfluid_gpu_impl, dimension>::density(float_type value)
{
    // set particle density
    density_ = value;
    LOG("particle density: " << density_);

    // compute periodic box length
    box_ = powf(npart / density_, 1. / dimension);
    LOG("periodic simulation box length: " << box_);

    // copy periodic box length to device symbol
    ljfluid_gpu_impl<dimension>::box(box_);
}

template <template <int> class ljfluid_gpu_impl, int dimension>
void ljfluid_gpu<ljfluid_gpu_impl, dimension>::box(float_type value)
{
    // set periodic box length
    box_ = value;
    LOG("periodic simulation box length: " << box_);

    // compute particle density
    density_ = npart / powf(box_, dimension);
    LOG("particle density: " << density_);

    // copy periodic box length to device symbol
    ljfluid_gpu_impl<dimension>::box(box_);
}

template <template <int> class ljfluid_gpu_impl, int dimension>
void ljfluid_gpu<ljfluid_gpu_impl, dimension>::timestep(float_type value)
{
    // set simulation timestep
    timestep_ = value;
    LOG("simulation timestep: " << timestep_);
    // copy simulation timestep to device symbol
    ljfluid_gpu_impl<dimension>::timestep(timestep_);
}

/**
 * seed random number generator
 */
template <template <int> class ljfluid_gpu_impl, int dimension>
void ljfluid_gpu<ljfluid_gpu_impl, dimension>::rng(unsigned int seed)
{
    LOG("random number generator seed: " << seed);
    try {
	rng_.set(seed, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const&) {
	throw exception("failed to seed random number generator");
    }
}

/**
 * restore random number generator from state
 */
template <template <int> class ljfluid_gpu_impl, int dimension>
void ljfluid_gpu<ljfluid_gpu_impl, dimension>::rng(rand48::state_type const& state)
{
    try {
	rng_.restore(state);
    }
    catch (cuda::error const&) {
	throw exception("failed to restore random number generator state");
    }
}

template <template <int> class ljfluid_gpu_impl, int dimension>
void ljfluid_gpu<ljfluid_gpu_impl, dimension>::attrs(H5::Group const& param) const
{
    H5xx::group node(param.openGroup("mdsim"));
    node["blocks"] = dim_.blocks_per_grid();
    node["threads"] = dim_.threads_per_block();

    // implementation-dependent attributes
    ljfluid_gpu_impl<dimension>::attrs(param);
}

// implementations for use with ljfluid class template
template <int dimension>
class ljfluid_gpu_square :
    public ljfluid_gpu<ljfluid_gpu_impl_square, dimension> {};

template <int dimension>
class ljfluid_gpu_cell :
    public ljfluid_gpu<ljfluid_gpu_impl_cell, dimension> {};

template <int dimension>
class ljfluid_gpu_neighbour :
    public ljfluid_gpu<ljfluid_gpu_impl_neighbour, dimension> {};

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_GPU_HPP */
