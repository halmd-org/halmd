/* Block correlations
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

#ifndef MDSIM_CORRELATION_HPP
#define MDSIM_CORRELATION_HPP

#include <H5Cpp.h>
#include <algorithm>
#include <boost/array.hpp>
// requires patch from http://svn.boost.org/trac/boost/ticket/1852
#include <boost/circular_buffer.hpp>
#include <boost/multi_array.hpp>
#include <boost/variant.hpp>
#include <cmath>
#include <string>
#include <vector>
#include "H5param.hpp"
#include "config.hpp"
#include "tcf.hpp"

namespace mdsim {

/**
 * Phase space sample
 */
struct phase_space_point
{
    typedef std::vector<hvector> vector_type;
    typedef hvector::value_type value_type;
    typedef std::vector<std::pair<hvector, hvector> > density_vector_type;

    phase_space_point(vector_type const& r, vector_type const& v, std::vector<value_type> const& q) : r(r), v(v), rho(q.size(), std::pair<hvector, hvector>(0, 0))
    {
	// spatial Fourier transformation
 	for (size_t i = 0; i < r.size(); ++i) {
	    for (unsigned int j = 0; j < q.size(); ++j) {
		// compute averages to maintain accuracy with single precision floating-point
		rho[j].first += (cos(r[i] * q[j]) - rho[j].first) / (i + 1);
		rho[j].second += (sin(r[i] * q[j]) - rho[j].second) / (i + 1);
	    }
	}
	// normalize Fourier transformed density with N^(-1/2)
	const value_type n = std::sqrt(r.size());
	for (unsigned int j = 0; j < q.size(); ++j) {
	    // multiply averages with N^(+1/2)
	    rho[j].first *= n;
	    rho[j].second *= n;
	}
    }

    /** particle positions */
    vector_type r;
    /** particle velocities */
    vector_type v;
    /** spatially Fourier transformed density for given q-values */
    density_vector_type rho;
};

/**
 * Block correlations
 */
class correlation
{
public:
    /** phase space sample type */
    typedef phase_space_point sample_type;
    /** phase space samples block type */
    typedef boost::circular_buffer<phase_space_point> block_type;
    /** correlation function variant type */
    typedef boost::make_variant_over<tcf_types>::type tcf_type;

public:
    /** initialize with all correlation function types */
    correlation();
    /** set total number of simulation steps */
    void steps(uint64_t const& value, double const& timestep);
    /** set total simulation time */
    void time(double const& value, double const& timestep);
    /** set sample rate for lowest block level */
    void sample_rate(unsigned int const& sample_rate);
    /** set block size */
    void block_size(unsigned int const& value);
    /** set maximum number of samples per block */
    void max_samples(uint64_t const& value);
    /** set q-vectors for spatial Fourier transformation */
    void q_values(unsigned int const& n, double const& box);
 
    /** returns total number of simulation steps */
    uint64_t const& steps() const { return m_steps; }
    /** returns total simulation time */
    double const& time() const { return m_time; }
    /** returns sample rate for lowest block level */
    unsigned int const& sample_rate() const { return m_sample_rate; }
    /** returns block size */
    unsigned int const& block_size() const { return m_block_size; }
    /** returns block shift */
    unsigned int const& block_shift() const { return m_block_shift; }
    /** returns block count */
    unsigned int const& block_count() const { return m_block_count; }
    /** returns maximum number of samples per block */
    uint64_t const& max_samples() const { return m_max_samples; }
    /** returns number of wave vectors */
    unsigned int q_values() { return m_q_vector.size(); }

    /** create HDF5 correlations output file */
    void open(std::string const& filename);
    /** close HDF5 file */
    void close();
    /** returns HDF5 parameter group */
    H5param attrs();
    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;
    /** check if sample is acquired for given simulation step */
    bool sample(uint64_t const& step) const;
    /** sample time correlation functions */
    void sample(std::vector<hvector> const& r, std::vector<hvector> const& v, uint64_t const& step, bool& flush);
    /** write correlation function results to HDF5 file */
    void flush();

private:
    /** apply correlation functions to block samples */
    void autocorrelate_block(unsigned int n);

private:
    /** phase space sample blocks */
    std::vector<block_type> m_block;
    /** phase sample frequencies for block levels */
    std::vector<uint64_t> m_block_freq;
    /** correlation sample counts for block levels */
    std::vector<uint64_t> m_block_samples;

    /** simulation timestep */
    double m_timestep;
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
    boost::multi_array<double, 2> m_block_time;
    /** maximum number of correlation samples per block */
    uint64_t m_max_samples;
    /** q-values for spatial Fourier transformation */
    std::vector<double> m_q_vector;

    /** correlation functions and results */
    std::vector<tcf_type> m_tcf;
    /** HDF5 output file */
    H5::H5File m_file;
};

} // namespace mdsim

#endif /* ! MDSIM_AUTOCORRELATION_HPP */
