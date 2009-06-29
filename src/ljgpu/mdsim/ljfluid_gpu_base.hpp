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

#include <boost/assign/list_of.hpp>
#include <cuda_wrapper.hpp>
#include <ljgpu/algorithm/radix_sort.hpp>
#include <ljgpu/algorithm/reduce.hpp>
#include <ljgpu/math/gpu/dsfloat.cuh>
#include <ljgpu/mdsim/gpu/boltzmann.hpp>
#include <ljgpu/mdsim/gpu/lattice.hpp>
#include <ljgpu/mdsim/gpu/ljfluid_cell.hpp>
#include <ljgpu/mdsim/gpu/ljfluid_nbr.hpp>
#include <ljgpu/mdsim/gpu/ljfluid_square.hpp>
#include <ljgpu/mdsim/ljfluid_base.hpp>
#include <ljgpu/mdsim/virial.hpp>
#include <ljgpu/rng/rand48.hpp>
#include <ljgpu/sample/H5param.hpp>

namespace ljgpu
{

template <typename ljfluid_impl, int dimension>
class ljfluid_gpu_base : public ljfluid_base<ljfluid_impl, dimension>
{
public:
    typedef ljfluid_base<ljfluid_impl, dimension> _Base;
    typedef gpu::ljfluid<ljfluid_impl, dimension> _gpu;
    typedef typename _Base::float_type float_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::traits_type::gpu_vector_type gpu_vector_type;
    typedef typename _Base::host_sample_type host_sample_type;
    typedef typename _Base::energy_sample_type energy_sample_type;
    typedef typename _Base::virial_tensor virial_tensor;
    typedef typename _Base::traits_type::gpu_sample_type gpu_sample_type;

public:
    /** set number of particles in system */
    template <typename T>
    void particles(T const& value);
    /** set potential well depths */
    void epsilon(boost::array<float, 3> const& value);
    /** set collision diameters */
    void sigma(boost::array<float, 3> const& value);
    /** set number of CUDA execution threads */
    void threads(unsigned int value);
    /* set particle density */
    void density(float_type value);
    /* set periodic box length */
    void box(float_type value);
    /** set simulation timestep */
    void timestep(double value);
    /** set potential cutoff radius */
    void cutoff_radius(float_type value);
    /** set potential smoothing function scale parameter */
    void potential_smoothing(float_type value);

    /** seed random number generator */
    void rng(unsigned int seed);
    /** restore random number generator from state */
    void rng(rand48::state_type const& state);

    /** returns number of particles */
    unsigned int particles() const { return npart; }
    /** returns particle density */
    float_type density() const { return density_; }
    /** returns periodic box length */
    float_type box() const { return box_; }
    /** returns simulation timestep */
    double timestep() const { return timestep_; }
    /** returns potential cutoff radius */
    float_type cutoff_radius() const { return r_cut_sigma; }

    /** get number of CUDA execution blocks */
    unsigned int blocks() const { return dim_.blocks_per_grid(); }
    /** get number of CUDA execution threads */
    unsigned int threads() const { return dim_.threads_per_block(); }

    /** write parameters to HDF5 parameter group */
    void param(H5param& param) const;

protected:
    /** place particles on an fcc lattice */
    void lattice(cuda::vector<float4>& g_r);
    /** randomly permute particle coordinates */
    void random_permute(cuda::vector<float4>& g_r);
    /** assign ascending particle numbers */
    void init_tags(cuda::vector<float4>& g_r, cuda::vector<unsigned int>& g_tag);
    /** rescale particle velocities */
    void rescale_velocities(cuda::vector<gpu_vector_type>& g_v,
			double coeff,
			cuda::config const& dim,
			cuda::stream& stream);
    /** generate Maxwell-Boltzmann distributed velocities */
    void boltzmann(cuda::vector<gpu_vector_type>& g_v, double temp, cuda::stream& stream);
    /** restore system state from phase space sample */
    void state(host_sample_type& sample, float_type box,
	       cuda::host::vector<float4>& h_r,
	       cuda::host::vector<gpu_vector_type>& h_v);
    /** update Lennard-Jones forces */
    void update_forces(cuda::vector<float4>& r,
		       cuda::vector<gpu_vector_type>& v,
		       cuda::vector<gpu_vector_type>& f,
		       cuda::vector<float>& en,
		       cuda::vector<gpu_vector_type>& virial);

protected:
    using _Base::box_;
    using _Base::density_;
    using _Base::en_cut;
    using _Base::epsilon_;
    using _Base::m_times;
    using _Base::mixture_;
    using _Base::mpart;
    using _Base::npart;
    using _Base::potential_;
    using _Base::r_cut;
    using _Base::r_cut_sigma;
    using _Base::r_smooth;
    using _Base::rr_cut;
    using _Base::rri_smooth;
    using _Base::sigma2_;
    using _Base::timestep_;

    /** CUDA execution dimensions */
    cuda::config dim_;
    /** CUDA asynchronous execution */
    cuda::stream mutable stream_;
    /** CUDA timing */
    boost::array<cuda::event, 2> mutable event_;
    /** GPU random number generator */
    rand48 rng_;
    /** block sum of velocity */
    cuda::vector<gpu_vector_type> g_vcm;
    /** block sum of squared velocity */
    cuda::vector<dsfloat> g_vv;
    /** GPU radix sort */
    radix_sort<float4> radix_sort_;
    /** center of mass velocity */
    reduce<tag::sum, gpu_vector_type, vector_type> mutable reduce_velocity;
    /** squared velocity sum */
    reduce<tag::sum_of_squares, dsfloat, double> mutable reduce_squared_velocity;
    /** potential energy sum */
    reduce<tag::sum, dsfloat, double> mutable reduce_en;
    /** virial equation sum */
    virial_sum<ljgpu::cu::vector<dsfloat, virial_tensor::static_size>, virial_tensor> mutable reduce_virial;
};

template <typename ljfluid_impl, int dimension>
template <typename T>
void ljfluid_gpu_base<ljfluid_impl, dimension>::particles(T const& value)
{
    _Base::particles(value);

    try {
	cuda::copy(npart, _gpu::npart);
	cuda::copy(mpart, _gpu::mpart);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to copy particle number to device symbol");
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::epsilon(boost::array<float, 3> const& value)
{
    _Base::epsilon(value);

    try {
	cuda::copy(epsilon_, _gpu::epsilon);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to copy collision diameters to device symbol");
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::sigma(boost::array<float, 3> const& value)
{
    _Base::sigma(value);

    try {
	cuda::copy(sigma2_, _gpu::sigma2);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to copy squared collision diameters to device symbol");
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::cutoff_radius(float_type value)
{
    _Base::cutoff_radius(value);

    try {
	cuda::copy(r_cut, _gpu::r_cut);
	cuda::copy(rr_cut, _gpu::rr_cut);
	cuda::copy(en_cut, _gpu::en_cut);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw potential_energy_divergence();
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::potential_smoothing(float_type value)
{
    _Base::potential_smoothing(value);

    try {
	cuda::copy(rri_smooth, _gpu::rri_smooth);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to copy potential smoothing function scale symbol");
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::threads(unsigned int value)
{
    // query CUDA device properties
    cuda::device::properties prop;
    try {
	prop = cuda::device::properties(cuda::device::get());
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
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

    try {
	radix_sort_.resize(npart, dim_.threads_per_block());
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to allocate global device memory for radix sort");
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::density(float_type value)
{
    _Base::density(value);

    try {
	cuda::copy(box_, _gpu::box);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to copy periodic box length to device symbol");
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::box(float_type value)
{
    _Base::box(value);

    try {
	cuda::copy(box_, _gpu::box);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to copy periodic box length to device symbol");
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::timestep(double value)
{
    _Base::timestep(value);

    try {
	cuda::copy(static_cast<float_type>(timestep_), _gpu::timestep);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to copy simulation timestep to device symbol");
    }
}

/**
 * seed random number generator
 */
template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::rng(unsigned int seed)
{
    typedef gpu::boltzmann<dimension> _gpu;

    try {
	rng_.resize(cuda::config(_gpu::BLOCKS, _gpu::THREADS));
#ifdef USE_VERLET_DSFUN
	g_vcm.resize(2 * _gpu::BLOCKS);
#else
	g_vcm.resize(_gpu::BLOCKS);
#endif
	g_vv.resize(_gpu::BLOCKS);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to change random number generator dimensions");
    }

    LOG("random number generator seed: " << seed);

    try {
	rng_.set(seed, stream_);
	stream_.synchronize();
	rng_.init_symbols(_gpu::rand48::a, _gpu::rand48::c, _gpu::rand48::state);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to seed random number generator");
    }
}

/**
 * restore random number generator from state
 */
template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::rng(rand48::state_type const& state)
{
    try {
	rng_.restore(state);
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to restore random number generator state");
    }
}

/**
 * place particles on an fcc lattice
 */
template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::lattice(cuda::vector<float4>& g_r)
{
    LOG("placing particles on face-centered cubic (fcc) lattice");

    // particles per 2- or 3-dimensional unit cell
    unsigned int const m = 2 * (dimension - 1);
    // lower boundary for number of particles per lattice dimension
    unsigned int n = static_cast<unsigned int>(std::pow(npart / m, 1.f / dimension));
    // lower boundary for total number of lattice sites
    unsigned int N = m * static_cast<unsigned int>(pow(n, static_cast<unsigned int>(dimension)));

    if (N < npart) {
	n += 1;
	N = m * static_cast<unsigned int>(pow(n, dimension));
    }
    if (N > npart) {
	LOG_WARNING("lattice not fully occupied (" << N << " sites)");
    }

    // minimum distance in 2- or 3-dimensional fcc lattice
    LOG("minimum lattice distance: " << (box_ / n) / std::sqrt(2.f));

    try {
	event_[0].record(stream_);
	cuda::configure(dim_.grid, dim_.block, stream_);
	gpu::lattice<dimension>::fcc(g_r, n, box_);
	event_[1].record(stream_);
	event_[1].synchronize();
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to compute particle lattice positions on GPU");
    }

    m_times["lattice"] += event_[1] - event_[0];
}

/**
 * randomly assign particle types in a binary mixture
 */
template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::random_permute(cuda::vector<float4>& g_r)
{
    LOG("randomly permuting particle coordinates");

    cuda::vector<unsigned int> g_sort_index(npart);
    g_sort_index.reserve(dim_.threads());

    try {
	rng_.get(g_sort_index, stream_);
	radix_sort_(g_sort_index, g_r, stream_);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to randomly permute particle coordinates on GPU");
    }
}

/**
 * assign ascending particle numbers
 */
template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::init_tags(cuda::vector<float4>& g_r,
					       cuda::vector<unsigned int>& g_tag)
{
    try {
	cuda::configure(dim_.grid, dim_.block, stream_);
	_gpu::init_tags(g_r, g_tag);
	stream_.synchronize();
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to initialise particle tags on GPU");
    }
}

/**
 * rescale particle velocities
 */
template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::rescale_velocities(cuda::vector<gpu_vector_type>& g_v,
						    double coeff,
						    cuda::config const& dim,
						    cuda::stream& stream)
{
    try {
	cuda::configure(dim.grid, dim.block, stream);
	_gpu::rescale_velocity(g_v, coeff);
	stream.synchronize();
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	throw exception("failed to rescale velocities on GPU");
    }
}

/**
 * generate Maxwell-Boltzmann distributed velocities
 */
template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::boltzmann(cuda::vector<gpu_vector_type>& g_v, double temp, cuda::stream& stream)
{
    typedef gpu::boltzmann<dimension> _gpu;

    //
    // The particle velocities need to fullfill two constraints:
    //
    //  1. center of mass velocity shall be zero
    //  2. temperature of the distribution shall equal exactly the given value
    //
    // The above order is chosen as shifting the center of mass velocity
    // means altering the first moment of the velocity distribution, which
    // in consequence affects the second moment, i.e. the temperature.
    //

    // generate Maxwell-Boltzmann distributed velocities and reduce velocity
    cuda::configure(_gpu::BLOCKS, _gpu::THREADS, stream);
    _gpu::gaussian(g_v, npart, dim_.threads(), temp, g_vcm);
    // set center of mass velocity to zero and reduce squared velocity
    cuda::configure(_gpu::BLOCKS, _gpu::THREADS, stream);
    _gpu::shift_velocity(g_v, npart, dim_.threads(), g_vcm, g_vv);
    // rescale velocities to accurate temperature
    cuda::configure(_gpu::BLOCKS, _gpu::THREADS, stream);
    _gpu::scale_velocity(g_v, npart, dim_.threads(), g_vv, temp);
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::state(host_sample_type& sample, float_type box,
					    cuda::host::vector<float4>& h_r,
					    cuda::host::vector<gpu_vector_type>& h_v)
{
    _Base::state(sample, box);

    for (size_t i = 0, n = 0; n < npart; n += mpart[i], ++i) {
	std::copy(sample[i].r->begin(), sample[i].r->end(), h_r.begin() + n);
	std::copy(sample[i].v->begin(), sample[i].v->end(), h_v.begin() + n);
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::update_forces(cuda::vector<float4>& r,
						   cuda::vector<gpu_vector_type>& v,
						   cuda::vector<gpu_vector_type>& f,
						   cuda::vector<float>& en,
						   cuda::vector<gpu_vector_type>& virial)
{
    // (CUDA kernel execution is configured in derived class)
    if (mixture_ == BINARY) {
	if (potential_ == C2POT) {
	    _gpu::template variant<BINARY, C2POT>::mdstep(r, v, f, en, virial);
	}
	else {
	    _gpu::template variant<BINARY, C0POT>::mdstep(r, v, f, en, virial);
	}
    }
    else {
	if (potential_ == C2POT) {
	    _gpu::template variant<UNARY, C2POT>::mdstep(r, v, f, en, virial);
	}
	else {
	    _gpu::template variant<UNARY, C0POT>::mdstep(r, v, f, en, virial);
	}
    }
}

template <typename ljfluid_impl, int dimension>
void ljfluid_gpu_base<ljfluid_impl, dimension>::param(H5param& param) const
{
    _Base::param(param);

    H5xx::group node(param["mdsim"]);
    node["blocks"] = blocks();
    node["threads"] = threads();
}

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_LJFLUID_GPU_BASE_HPP */
