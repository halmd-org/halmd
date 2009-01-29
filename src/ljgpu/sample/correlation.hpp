/* Block correlations
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

#ifndef LJGPU_SAMPLE_CORRELATION_HPP
#define LJGPU_SAMPLE_CORRELATION_HPP

#include <H5Cpp.h>
// requires boost 1.37.0 or patch from http://svn.boost.org/trac/boost/ticket/1852
#include <boost/circular_buffer.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/variant.hpp>
#include <cmath>
#ifdef WITH_CUDA
# include <ljgpu/sample/tcf_gpu.hpp>
#endif
#include <ljgpu/sample/tcf_host.hpp>
#include <ljgpu/sample/tcf_visitor.hpp>
#include <ljgpu/sample/H5param.hpp>
#include <ljgpu/util/log.hpp>
#include <string>
#include <vector>

namespace ljgpu {

/**
 * Block correlations
 */
template <int dimension>
class correlation
{
public:
    typedef std::vector<tcf_host_sample<dimension> > host_sample_type;
    typedef boost::circular_buffer<host_sample_type> host_block_type;
#ifdef WITH_CUDA
    typedef std::vector<tcf_gpu_sample<dimension> > gpu_sample_type;
    typedef boost::circular_buffer<gpu_sample_type> gpu_block_type;
    typedef boost::variant<gpu_sample_type, host_sample_type> sample_variant;
    typedef boost::variant<gpu_block_type, host_block_type> block_variant;
    typedef typename boost::mpl::insert_range<tcf_gpu_types, typename boost::mpl::end<tcf_gpu_types>::type, tcf_host_types>::type tcf_types;
#else
    typedef boost::variant<host_sample_type> sample_variant;
    typedef boost::variant<host_block_type> block_variant;
    typedef tcf_host_types tcf_types;
#endif
    typedef boost::make_variant_over<tcf_types>::type tcf_variant;
    typedef std::vector<tcf_variant> tcf_vector;

    typedef tcf_sample<dimension> sample_type;
    typedef typename sample_type::vector_type vector_type;
    typedef typename sample_type::q_value_vector q_value_vector;
    typedef typename sample_type::q_vector_vector q_vector_vector;
    typedef boost::multi_array<double, 2> block_time_type;

public:
    /** set total number of simulation steps */
    void steps(uint64_t value, float timestep);
    /** set total simulation time */
    void time(double value, float timestep);
    /** set sample rate for lowest block level */
    void sample_rate(unsigned int sample_rate);
    /** set block size */
    void block_size(unsigned int value);
    /** set maximum number of samples per block */
    void max_samples(uint64_t value);
    /** set q-vectors for spatial Fourier transformation */
    void q_values(std::vector<float> const& values, float error, float box);
#ifdef WITH_CUDA
    /** add correlation functions for GPU */
    void add_gpu_correlation_functions(size_t types);
#endif
    /** add correlation functions for host */
    void add_host_correlation_functions(size_t types);

    /** returns total number of simulation steps */
    uint64_t steps() const { return m_steps; }
    /** returns total simulation time */
    double time() const { return m_time; }
    /** returns sample rate for lowest block level */
    unsigned int sample_rate() const { return m_sample_rate; }
    /** returns block size */
    unsigned int block_size() const { return m_block_size; }
    /** returns block shift */
    unsigned int block_shift() const { return m_block_shift; }
    /** returns block count */
    unsigned int block_count() const { return m_block_count; }
    /** returns maximum number of samples per block */
    uint64_t max_samples() const { return m_max_samples; }
    /** returns number of wave vectors */
    unsigned int q_values() { return m_q_vector.size(); }

    /** create HDF5 correlations output file */
    void open(std::string const& filename, size_t types);
    /** close HDF5 file */
    void close();
    /** returns HDF5 parameter group */
    operator H5param() { return m_file; }
    /** write parameters to HDF5 parameter group */
    void param(H5::Group const& param) const;
    /** check if sample is acquired for given simulation step */
    bool sample(int64_t step) const;
    /** sample time correlation functions */
    template <typename trajectory_sample_variant>
    void sample(trajectory_sample_variant const& sample_, uint64_t step, bool& flush);
    /** write correlation function results to HDF5 file */
    void flush();

private:
    /** compute lattice points in first octant on surface of 3-dimensional spheres */
    template <typename T>
    typename boost::enable_if<boost::is_same<vector<double, 3>, T>, void>::type
    find_q_vectors(std::vector<std::pair<int, int> > const& qq, int q_max, std::vector<std::vector<T> >& q);
    /** compute lattice points in first quadrant on surface of 2-dimensional spheres */
    template <typename T>
    typename boost::enable_if<boost::is_same<vector<double, 2>, T>, void>::type
    find_q_vectors(std::vector<std::pair<int, int> > const& qq, int q_max, std::vector<std::vector<T> >& q);
    /** apply correlation functions to block samples */
    void autocorrelate_block(unsigned int n);

private:
    /** phase space sample blocks */
    std::vector<block_variant> m_block;
    /** GPU or host sample */
    sample_variant m_sample;
    /** phase sample frequencies for block levels */
    std::vector<uint64_t> m_block_freq;
    /** correlation sample counts for block levels */
    std::vector<uint64_t> m_block_samples;

    /** simulation timestep */
    float m_timestep;
    /** sample rate for lowest block level */
    unsigned int m_sample_rate;
    /** total number of simulation steps */
    uint64_t m_steps;
    /** total simulation time */
    double m_time;
    /** block size */
    unsigned int m_block_size;
    /** block shift */
    unsigned int m_block_shift;
    /** block count */
    unsigned int m_block_count;
    /** block time intervals */
    block_time_type m_block_time;
    /** maximum number of correlation samples per block */
    uint64_t m_max_samples;
    /** q values for spatial Fourier transformation */
    q_value_vector m_q_value;
    /** q vectors for spatial Fourier transformation */
    q_vector_vector m_q_vector;
    /** relative deviation of averaging wave vectors */
    float m_q_error;

    /** correlation functions and results */
    tcf_vector m_tcf;
    /** HDF5 output file */
    H5::H5File m_file;
};

/**
 * sample time correlation functions
 */
template <int dimension>
template <typename trajectory_sample_variant>
void correlation<dimension>::sample(trajectory_sample_variant const& sample_, uint64_t step, bool& flush)
{
    // empty sample of GPU or host type
    sample_variant sample(m_sample);
    // copy phase space coordinates and compute spatial Fourier transformation 
    boost::apply_visitor(tcf_sample_phase_space(), sample, sample_);
    boost::apply_visitor(tcf_fourier_transform_sample(m_q_vector), sample);

    for (unsigned int i = 0; i < m_block_count; ++i) {
	if (m_block_samples[i] >= m_max_samples)
	    continue;
	if (step % m_block_freq[i])
	    continue;

	boost::apply_visitor(tcf_block_add_sample(), m_block[i], sample);

	if (boost::apply_visitor(tcf_block_is_full(), m_block[i])) {
	    autocorrelate_block(i);
	    if (i < 2) {
		// sample_ only full blocks in lowest levels to account for strong correlations
		boost::apply_visitor(tcf_block_clear(), m_block[i]);
	    }
	    m_block_samples[i]++;

	    if (m_max_samples == m_block_samples[i]) {
		LOG("finished sampling on block level " << i << " at step " << step);
		// trigger global write of partial results to disk
		flush = true;
	    }
	}
    }
}

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_CORRELATION_HPP */
