/* Time correlation functions
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

#ifndef LJGPU_SAMPLE_TCF_HPP
#define LJGPU_SAMPLE_TCF_HPP

#include <algorithm>
#include <boost/array.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/multi_array.hpp>
#include <boost/variant.hpp>
#include <cuda_wrapper.hpp>
#include <ljgpu/math/accum.hpp>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/mdsim/sample.hpp>
#include <ljgpu/mdsim/traits.hpp>
#include <ljgpu/sample/gpu/tcf.hpp>
#include <vector>

namespace ljgpu {

/**
 * Phase space sample for evaluating correlation functions
 */
template <int dimension>
struct tcf_gpu_sample
{
    typedef vector<double, dimension> vector_type;
    typedef std::vector<vector_type> sample_vector;
    typedef std::vector<double> q_value_vector;
    typedef std::vector<std::vector<vector_type> > q_vector_vector;

    /** real and imaginary components of Fourier transformed density rho(q) */
    typedef std::pair<double, double> density_vector_pair;
    /** vector of Fourier transformed densities for different q-values */
    typedef std::vector<std::vector<density_vector_pair> > density_vector_vector;

    typedef mdsim_traits<ljfluid_impl_gpu_base<dimension> > traits_type;
    typedef typename traits_type::gpu_vector_type gpu_vector_type;
    typedef cuda::vector<gpu_vector_type> gpu_sample_vector;
    typedef gpu::tcf<dimension> _gpu;

    /**
     * initialise phase space sample
     */
    void operator()(q_vector_vector const& q)
    {
	// CUDA execution dimensions for accumulator
	enum { BLOCKS = _gpu::BLOCKS };
	enum { THREADS = _gpu::THREADS };

	// q-vector iterators
	typename q_vector_vector::const_iterator q0;
	typename q_vector_vector::value_type::const_iterator q1;
	// Fourier-transformed density iterators
	typename density_vector_vector::iterator rho0;
	typename density_vector_vector::value_type::iterator rho1;
	// accumulator iterator
	dfloat* sum0;

	// allocate memory for Fourier-transformed densities
	size_t size = 0;
	rho.resize(q.size());
	for (q0 = q.begin(), rho0 = rho.begin(); q0 != q.end(); ++q0, ++rho0) {
	    rho0->assign(q0->size(), density_vector_pair(0, 0));
	    size += q0->size();
	}
	// allocate device and host memory for accumulators
	cuda::vector<dfloat> g_sum(2 * size * BLOCKS);
	cuda::host::vector<dfloat> h_sum(g_sum.size());

	// compute Fourier-transformed densities on GPU
	for (q0 = q.begin(), sum0 = g_sum.data(); q0 != q.end(); ++q0) {
	    for (q1 = q0->begin(); q1 != q0->end(); ++q1, sum0 += 2 * BLOCKS) {
		cuda::configure(BLOCKS, THREADS);
		_gpu::coherent_scattering_function(r, *q1, sum0, sum0 + BLOCKS, r.size());
	    }
	}
	// copy accumulator block results from GPU to host
	cuda::copy(g_sum, h_sum);
	// accumulate Fourier-transformed densities on host
	for (rho0 = rho.begin(), sum0 = h_sum.data(); rho0 != rho.end(); ++rho0) {
	    for (rho1 = rho0->begin(); rho1 != rho0->end(); ++rho1, sum0 += 2 * BLOCKS) {
		rho1->first = std::accumulate(sum0, sum0 + BLOCKS, 0.);
		rho1->second = std::accumulate(sum0 + BLOCKS, sum0 + 2 * BLOCKS, 0.);
	    }
	}
    }

    /** particle positions */
    gpu_sample_vector r;
    /** particle velocities */
    gpu_sample_vector v;
    /** spatially Fourier transformed density for given q-values */
    density_vector_vector rho;
};

/** correlation function result types */
typedef boost::multi_array<accumulator<double>, 2> tcf_unary_result_type;
typedef boost::multi_array<accumulator<double>, 3> tcf_binary_result_type;

template <template <int> class sample_type>
struct correlation_function;

template <>
struct correlation_function<tcf_gpu_sample> {};

/**
 * mean-square displacement
 */
template <template <int> class sample_type>
struct mean_square_displacement;

template <>
struct mean_square_displacement<tcf_gpu_sample> : correlation_function<tcf_gpu_sample>
{
    /** block sample results */
    tcf_unary_result_type result;
    /** HDF5 dataset */
    H5::DataSet dataset;

    /** device and host memory for accumulators */
    cuda::vector<unsigned int> g_count;
    cuda::host::vector<unsigned int> h_count;
    cuda::vector<dfloat> g_mean;
    cuda::host::vector<dfloat> h_mean;
    cuda::vector<dfloat> g_var;
    cuda::host::vector<dfloat> h_var;

    char const* name() const { return "MSD"; }

    /**
     * autocorrelate samples in block
     */
    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type::value_type sample_type;
	typedef typename sample_type::_gpu _gpu;
	typedef typename output_iterator::value_type accumulator_type;
	enum { BLOCKS = _gpu::BLOCKS };
	enum { THREADS = _gpu::THREADS };

	sample_iterator sample;
	unsigned int *count, *count0;
	dfloat *mean, *var;

	// allocate device and host memory for accumulators, if necessary
	g_count.resize((last.first - first.first) * BLOCKS);
	h_count.resize(g_count.size());
	g_mean.resize((last.first - first.first) * BLOCKS);
	h_mean.resize(g_mean.size());
	g_var.resize((last.first - first.first) * BLOCKS);
	h_var.resize(g_var.size());

	// compute mean square displacements on GPU
	for (sample = first.first, count = g_count.data(), mean = g_mean.data(), var = g_var.data(); sample != last.first; ++sample, count += BLOCKS, mean += BLOCKS, var += BLOCKS) {
	    cuda::configure(BLOCKS, THREADS);
	    _gpu::mean_square_displacement((*sample)->r, (*first.first)->r, count, mean, var, (*sample)->r.size());
	}
	// copy accumulator block results from GPU to host
	cuda::copy(g_count, h_count);
	cuda::copy(g_mean, h_mean);
	cuda::copy(g_var, h_var);
	// accumulate mean square displacements on host
	for (sample = first.first, count0 = h_count.data(), mean = h_mean.data(), var = h_var.data(); sample != last.first; ++sample, ++result, count0 += BLOCKS) {
	    for (count = count0; count != count0 + BLOCKS; ++count, ++mean, ++var) {
		*result += accumulator_type(*count, *mean, *var);
	    }
	}
    }
};

/**
 * mean-quartic displacement
 */
template <template <int> class sample_type>
struct mean_quartic_displacement;

template <>
struct mean_quartic_displacement<tcf_gpu_sample> : correlation_function<tcf_gpu_sample>
{
    /** block sample results */
    tcf_unary_result_type result;
    /** HDF5 dataset */
    H5::DataSet dataset;

    /** device and host memory for accumulators */
    cuda::vector<unsigned int> g_count;
    cuda::host::vector<unsigned int> h_count;
    cuda::vector<dfloat> g_mean;
    cuda::host::vector<dfloat> h_mean;
    cuda::vector<dfloat> g_var;
    cuda::host::vector<dfloat> h_var;

    char const* name() const { return "MQD"; }

    /**
     * autocorrelate samples in block
     */
    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type::value_type sample_type;
	typedef typename sample_type::_gpu _gpu;
	typedef typename output_iterator::value_type accumulator_type;
	enum { BLOCKS = _gpu::BLOCKS };
	enum { THREADS = _gpu::THREADS };

	sample_iterator sample;
	unsigned int *count, *count0;
	dfloat *mean, *var;

	// allocate device and host memory for accumulators, if necessary
	g_count.resize((last.first - first.first) * BLOCKS);
	h_count.resize(g_count.size());
	g_mean.resize((last.first - first.first) * BLOCKS);
	h_mean.resize(g_mean.size());
	g_var.resize((last.first - first.first) * BLOCKS);
	h_var.resize(g_var.size());

	// compute mean quartic displacements on GPU
	for (sample = first.first, count = g_count.data(), mean = g_mean.data(), var = g_var.data(); sample != last.first; ++sample, count += BLOCKS, mean += BLOCKS, var += BLOCKS) {
	    cuda::configure(BLOCKS, THREADS);
	    _gpu::mean_quartic_displacement((*sample)->r, (*first.first)->r, count, mean, var, (*sample)->r.size());
	}
	// copy accumulator block results from GPU to host
	cuda::copy(g_count, h_count);
	cuda::copy(g_mean, h_mean);
	cuda::copy(g_var, h_var);
	// accumulate mean quartic displacements on host
	for (sample = first.first, count0 = h_count.data(), mean = h_mean.data(), var = h_var.data(); sample != last.first; ++sample, ++result, count0 += BLOCKS) {
	    for (count = count0; count != count0 + BLOCKS; ++count, ++mean, ++var) {
		*result += accumulator_type(*count, *mean, *var);
	    }
	}
    }
};

/**
 * velocity autocorrelation
 */
template <template <int> class sample_type>
struct velocity_autocorrelation;

template <>
struct velocity_autocorrelation<tcf_gpu_sample> : correlation_function<tcf_gpu_sample>
{
    /** block sample results */
    tcf_unary_result_type result;
    /** HDF5 dataset */
    H5::DataSet dataset;

    /** device and host memory for accumulators */
    cuda::vector<unsigned int> g_count;
    cuda::host::vector<unsigned int> h_count;
    cuda::vector<dfloat> g_mean;
    cuda::host::vector<dfloat> h_mean;
    cuda::vector<dfloat> g_var;
    cuda::host::vector<dfloat> h_var;

    char const* name() const { return "VAC"; }

    /**
     * autocorrelate samples in block
     */
    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type::value_type sample_type;
	typedef typename sample_type::_gpu _gpu;
	typedef typename output_iterator::value_type accumulator_type;
	enum { BLOCKS = _gpu::BLOCKS };
	enum { THREADS = _gpu::THREADS };

	sample_iterator sample;
	unsigned int *count, *count0;
	dfloat *mean, *var;

	// allocate device and host memory for accumulators, if necessary
	g_count.resize((last.first - first.first) * BLOCKS);
	h_count.resize(g_count.size());
	g_mean.resize((last.first - first.first) * BLOCKS);
	h_mean.resize(g_mean.size());
	g_var.resize((last.first - first.first) * BLOCKS);
	h_var.resize(g_var.size());

	// compute velocity autocorrelations on GPU
	for (sample = first.first, count = g_count.data(), mean = g_mean.data(), var = g_var.data(); sample != last.first; ++sample, count += BLOCKS, mean += BLOCKS, var += BLOCKS) {
	    cuda::configure(BLOCKS, THREADS);
	    _gpu::velocity_autocorrelation((*sample)->v, (*first.first)->v, count, mean, var, (*sample)->v.size());
	}
	// copy accumulator block results from GPU to host
	cuda::copy(g_count, h_count);
	cuda::copy(g_mean, h_mean);
	cuda::copy(g_var, h_var);
	// accumulate velocity autocorrelations on host
	for (sample = first.first, count0 = h_count.data(), mean = h_mean.data(), var = h_var.data(); sample != last.first; ++sample, ++result, count0 += BLOCKS) {
	    for (count = count0; count != count0 + BLOCKS; ++count, ++mean, ++var) {
		*result += accumulator_type(*count, *mean, *var);
	    }
	}
    }
};

/**
 * intermediate scattering function
 */
template <template <int> class sample_type>
struct intermediate_scattering_function;

template <>
struct intermediate_scattering_function<tcf_gpu_sample> : correlation_function<tcf_gpu_sample>
{
    /** block sample results */
    tcf_binary_result_type result;
    /** HDF5 dataset */
    H5::DataSet dataset;

    char const* name() const { return "ISF"; }

    /**
     * autocorrelate samples in block
     */
    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type::value_type sample_type;
	typedef typename sample_type::density_vector_vector density_vector_vector;
	typedef typename output_iterator::value_type result_vector;

	sample_iterator sample;
	typename density_vector_vector::const_iterator rho0, rho0_0;
	typename density_vector_vector::value_type::const_iterator rho1, rho1_0;
	typename result_vector::iterator result0;

	// accumulate intermediate scattering functions on host
	for (sample = first.first; sample != last.first; ++sample, ++result) {
	    for (rho0 = (*sample)->rho.begin(), rho0_0 = (*first.first)->rho.begin(), result0 = result->begin(); rho0 != (*sample)->rho.end(); ++rho0, ++rho0_0, ++result0) {
		for (rho1 = rho0->begin(), rho1_0 = rho0_0->begin(); rho1 != rho0->end(); ++rho1, ++rho1_0) {
		    *result0 += (rho1->first * rho1_0->first + rho1->second * rho1_0->second) / (*sample)->r.size();
		}
	    }
	}
    }
};

/**
 * self-intermediate scattering function
 */
template <template <int> class sample_type>
struct self_intermediate_scattering_function;

template <>
struct self_intermediate_scattering_function<tcf_gpu_sample> : correlation_function<tcf_gpu_sample>
{
    /** block sample results */
    tcf_binary_result_type result;
    /** HDF5 dataset */
    H5::DataSet dataset;

    /** device and host memory for accumulators */
    cuda::vector<dfloat> g_sum;
    cuda::host::vector<dfloat> h_sum;

    char const* name() const { return "SISF"; }

    /**
     * autocorrelate samples in block
     */
    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type::value_type sample_type;
	typedef typename sample_type::_gpu _gpu;
	typedef typename sample_type::q_vector_vector q_vector_vector;
	typedef typename output_iterator::value_type result_vector;
	enum { BLOCKS = _gpu::BLOCKS };
	enum { THREADS = _gpu::THREADS };

	sample_iterator sample;
	typename q_vector_vector::const_iterator q0;
	typename q_vector_vector::value_type::const_iterator q1;
	typename result_vector::iterator result0;
	dfloat* sum;

	// allocate device and host memory for accumulators, if necessary
	size_t size = 0;
	for (q0 = first.second; q0 != last.second; ++q0) {
	    size += q0->size();
	}
	g_sum.resize((last.first - first.first) * size * BLOCKS);
	h_sum.resize(g_sum.size());

	// compute self-intermediate scattering functions on GPU
	for (sample = first.first, sum = g_sum.data(); sample != last.first; ++sample) {
	    for (q0 = first.second; q0 != last.second; ++q0) {
		for (q1 = q0->begin(); q1 != q0->end(); ++q1, sum += BLOCKS) {
		    cuda::configure(BLOCKS, THREADS);
		    _gpu::incoherent_scattering_function((*sample)->r, (*first.first)->r, *q1, sum, (*sample)->r.size());
		}
	    }
	}
	// copy accumulator block results from GPU to host
	cuda::copy(g_sum, h_sum);
	// accumulate self-intermediate scattering functions on host
	for (sample = first.first, sum = h_sum.data(); sample != last.first; ++sample, ++result) {
	    for (q0 = first.second, result0 = result->begin(); q0 != last.second; ++q0, ++result0) {
		for (q1 = q0->begin(); q1 != q0->end(); ++q1, sum += BLOCKS) {
		    *result0 += std::accumulate(sum, sum + BLOCKS, 0.) / (*sample)->r.size();
		}
	    }
	}
    }
};

/** correlation function types */
typedef boost::mpl::vector<mean_square_displacement<tcf_gpu_sample> > _tcf_types_0;
typedef boost::mpl::push_back<_tcf_types_0, mean_quartic_displacement<tcf_gpu_sample> >::type _tcf_types_1;
typedef boost::mpl::push_back<_tcf_types_1, velocity_autocorrelation<tcf_gpu_sample> >::type _tcf_types_2;
typedef boost::mpl::push_back<_tcf_types_2, intermediate_scattering_function<tcf_gpu_sample> >::type _tcf_types_3;
typedef boost::mpl::push_back<_tcf_types_3, self_intermediate_scattering_function<tcf_gpu_sample> >::type tcf_types;

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_TCF_HPP */
