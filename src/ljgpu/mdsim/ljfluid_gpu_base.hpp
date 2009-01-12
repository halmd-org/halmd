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

#ifndef LJGPU_MDSIM_LJFLUID_GPU_BASE_HPP
#define LJGPU_MDSIM_LJFLUID_GPU_BASE_HPP

#include <cuda_wrapper.hpp>
#include <ljgpu/mdsim/ljfluid_base.hpp>
#include <ljgpu/mdsim/gpu/ljfluid_cell.hpp>
#include <ljgpu/mdsim/gpu/ljfluid_nbr.hpp>
#include <ljgpu/mdsim/gpu/ljfluid_square.hpp>
#include <ljgpu/rng/rand48.hpp>

namespace ljgpu
{

template <typename ljfluid_impl>
class ljfluid_gpu_base : public ljfluid_base<ljfluid_impl>
{
public:
    typedef ljfluid_base<ljfluid_impl> _Base;
    typedef gpu::ljfluid<ljfluid_impl> _gpu;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::traits_type::gpu_vector_type gpu_vector_type;
    typedef typename _Base::sample_type sample_type;
    typedef typename sample_type::sample_visitor sample_visitor;
    enum { dimension = _Base::dimension };

public:
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
    /** set potential cutoff radius */
    void cutoff_radius(float_type value);
#ifdef USE_POTENTIAL_SMOOTHING
    /** set potential smoothing function scale parameter */
    void potential_smoothing(float_type value);
#endif

    /** seed random number generator */
    void rng(unsigned int seed);
    /** restore random number generator from state */
    void rng(rand48::state_type const& state);

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
    /** returns potential cutoff radius */
    float_type cutoff_radius() const { return r_cut; }

    /** write parameters to HDF5 parameter group */
    void param(H5::Group const& param) const;

protected:
    /** CUDA execution dimensions */
    cuda::config dim_;
    /** CUDA asynchronous execution */
    cuda::stream stream_;
    /** GPU random number generator */
    rand48 rng_;

    using _Base::npart;
    using _Base::density_;
    using _Base::box_;
    using _Base::timestep_;
    using _Base::r_cut;
    using _Base::rr_cut;
    using _Base::en_cut;
#ifdef USE_POTENTIAL_SMOOTHING
    using _Base::r_smooth;
    using _Base::rri_smooth;
#endif

    using _Base::m_sample;
    using _Base::m_times;
};

template <typename ljfluid_impl>
void ljfluid_gpu_base<ljfluid_impl>::particles(unsigned int value)
{
    _Base::particles(value);

    // allocate swappable host memory for trajectory sample
    try {
	m_sample.r.resize(npart);
	m_sample.v.resize(npart);
    }
    catch (cuda::error const&) {
	throw exception("failed to allocate swappable host memory for trajectory sample");
    }

    try {
	cuda::copy(npart, _gpu::npart);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy particle number to device symbol");
    }
}

template <typename ljfluid_impl>
void ljfluid_gpu_base<ljfluid_impl>::cutoff_radius(float_type value)
{
    _Base::cutoff_radius(value);

    try {
	cuda::copy(r_cut, _gpu::r_cut);
	cuda::copy(rr_cut, _gpu::rr_cut);
	cuda::copy(en_cut, _gpu::en_cut);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy potential cutoff symbols");
    }
}

#ifdef USE_POTENTIAL_SMOOTHING
template <typename ljfluid_impl>
void ljfluid_gpu_base<ljfluid_impl>::potential_smoothing(float_type value)
{
    _Base::potential_smoothing(value);

    try {
	cuda::copy(rri_smooth, _gpu::rri_smooth);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy potential smoothing function scale symbol");
    }
}
#endif /* USE_POTENTIAL_SMOOTHING */


template <typename ljfluid_impl>
void ljfluid_gpu_base<ljfluid_impl>::threads(unsigned int value)
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
}

template <typename ljfluid_impl>
void ljfluid_gpu_base<ljfluid_impl>::density(float_type value)
{
    _Base::density(value);

    try {
	cuda::copy(box_, _gpu::box);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy periodic box length to device symbol");
    }
}

template <typename ljfluid_impl>
void ljfluid_gpu_base<ljfluid_impl>::box(float_type value)
{
    _Base::box(value);

    try {
	cuda::copy(box_, _gpu::box);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy periodic box length to device symbol");
    }
}

template <typename ljfluid_impl>
void ljfluid_gpu_base<ljfluid_impl>::timestep(float_type value)
{
    _Base::timestep(value);

    try {
	cuda::copy(timestep_, _gpu::timestep);
    }
    catch (cuda::error const&) {
	throw exception("failed to copy simulation timestep to device symbol");
    }
}

/**
 * seed random number generator
 */
template <typename ljfluid_impl>
void ljfluid_gpu_base<ljfluid_impl>::rng(unsigned int seed)
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
template <typename ljfluid_impl>
void ljfluid_gpu_base<ljfluid_impl>::rng(rand48::state_type const& state)
{
    try {
	rng_.restore(state);
    }
    catch (cuda::error const&) {
	throw exception("failed to restore random number generator state");
    }
}

template <typename ljfluid_impl>
void ljfluid_gpu_base<ljfluid_impl>::param(H5::Group const& param) const
{
    _Base::param(param);

    H5xx::group node(param.openGroup("mdsim"));
    node["blocks"] = blocks();
    node["threads"] = threads();
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_GPU_BASE_HPP */
