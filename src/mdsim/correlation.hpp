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
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/variant.hpp>
#include <cmath>
#include <string>
#include <unistd.h>
#include <vector>
#include "H5param.hpp"
#include "H5xx.hpp"
#include "accumulator.hpp"
#include "exception.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "tcf.hpp"

#define foreach BOOST_FOREACH

namespace mdsim {

/**
 * Block correlations
 */
template <unsigned dimension, typename T>
class correlation
{
public:
    /** phase space sample type */
    typedef phase_space_point<dimension, T> sample_type;
    /** phase space samples block type */
    typedef boost::circular_buffer<phase_space_point<dimension, T> > block_type;

    /** correlation function type */
    typedef typename boost::make_variant_over<tcf_types>::type tcf_type;
    /** correlation function results type */
    typedef boost::multi_array<accumulator<double>, 2> tcf_result_type;
    /** correlation function and results pair type */
    typedef std::pair<tcf_type, tcf_result_type> tcf_pair;

    /** binary correlation function type */
    typedef typename boost::make_variant_over<qtcf_types>::type qtcf_type;
    /** binary correlation function results type */
    typedef boost::multi_array<accumulator<double>, 3> qtcf_result_type;
    /** binary correlation function and results pair type */
    typedef std::pair<qtcf_type, qtcf_result_type> qtcf_pair;

public:
    /** set total number of simulation steps */
    void steps(uint64_t const& value, double const& timestep);
    /** set total simulation time */
    void time(double const& value, double const& timestep);
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
    /** returns sample timestep */
    double const& timestep() const { return m_timestep; }
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
    void sample(std::vector<T> const& r, std::vector<T> const& v, uint64_t const& step);
    /** write correlation function results to HDF5 file */
    void write();

private:
    /** apply correlation functions to block samples */
    void autocorrelate_block(unsigned int n);
    /** returns simulation time belonging to block sample at given block offset */
    double block_sample_time(unsigned int block, unsigned int offset) const;

private:
    /** phase space sample blocks */
    std::vector<block_type> m_block;
    /** phase sample frequencies for block levels */
    std::vector<uint64_t> m_block_freq;
    /** correlation sample counts for block levels */
    std::vector<uint64_t> m_block_samples;

    /** simulation timestep */
    double m_timestep;
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
    /** maximum number of correlation samples per block */
    uint64_t m_max_samples;
    /** q-values for spatial Fourier transformation */
    std::vector<double> m_q_vector;

    /** correlation functions and results */
    boost::array<tcf_pair, 3> m_tcf;
    /** binary correlation functions and results */
    boost::array<qtcf_pair, 2> m_qtcf;

    /** HDF5 output file */
    H5::H5File m_file;
};

/**
 * set total number of simulation steps
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::steps(uint64_t const& value, double const& timestep)
{
    // set total number of simulation steps
    m_steps = value;
    LOG("total number of simulation steps: " << m_steps);
    // set simulation timestep
    m_timestep = timestep;
    LOG("sample timestep: " << m_timestep);
    // derive total simulation time
    m_time = value * m_timestep;
    LOG("total simulation time: " << m_time);
}

/**
 * set total simulation time
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::time(double const& value, double const& timestep)
{
    // set total simulation time
    m_time = value;
    LOG("total simulation time: " << m_time);
    // set simulation timestep
    m_timestep = timestep;
    LOG("sample timestep: " << m_timestep);
    // derive total number of simulation steps
    m_steps = roundf(m_time / m_timestep);
    LOG("total number of simulation steps: " << m_steps);
}

/**
 * set block size
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::block_size(unsigned int const& value)
{
    // set block size
    m_block_size = value;
    LOG("block size: " << m_block_size);

    // derive block shift from block size
    m_block_shift = std::floor(std::sqrt(m_block_size));
    LOG("block shift: " << m_block_shift);
    if (m_block_shift < 2) {
	throw exception("computed block shift is less than 2, larger block size required");
    }

    // derive block count from block size and block shift
    m_block_count = 0;
    for (unsigned int n = m_block_size; n <= m_steps; n *= m_block_size) {
	m_block_count++;
	if ((n * m_block_shift) > m_steps) {
	    break;
	}
	m_block_count ++;
    }
    LOG("block count: " << m_block_count);
    if (!m_block_count) {
	throw exception("computed block count is zero, more simulations steps required");
    }

    // allocate phase space sample blocks
    try {
	m_block.resize(m_block_count, block_type(m_block_size));
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate phase space sample blocks");
    }
    m_block_samples.resize(m_block_count, 0);

    // calculate phase sample frequencies
    for (unsigned int i = 0, n = 1, m = m_block_shift; i < m_block_count; ++i) {
	if (i % 2) {
	    m_block_freq.push_back(m);
	    m = m * m_block_size;
	}
	else {
	    m_block_freq.push_back(n);
	    n = n * m_block_size;
	}
    }

    // setup correlation functions
    m_tcf[0].first = mean_square_displacement();
    m_tcf[1].first = mean_quartic_displacement();
    m_tcf[2].first = velocity_autocorrelation();

    // allocate correlation functions results
    try {
	foreach (tcf_pair& tcf, m_tcf) {
	    tcf.second.resize(boost::extents[m_block_count][m_block_size - 1]);
	}
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate correlation functions results");
    }
}

/**
 * set maximum number of samples per block
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::max_samples(uint64_t const& value)
{
    m_max_samples = value;
    LOG("maximum number of samples per block: " << m_max_samples);
    if (m_max_samples < m_block_size) {
	throw exception("maximum number of samples must not be smaller than block size");
    }
}

/**
 * set q-vectors for spatial Fourier transformation
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::q_values(unsigned int const& n, double const& box)
{
    // integer multiples of q-value corresponding to periodic box length
    for (unsigned int k = 1; k <= n; ++k) {
	m_q_vector.push_back(k * 2 * M_PI / box);
    }

    // setup binary correlation functions
    m_qtcf[0].first = intermediate_scattering_function();
    m_qtcf[1].first = self_intermediate_scattering_function();

    // allocate binary correlation functions results
    try {
	foreach (qtcf_pair& qtcf, m_qtcf) {
	    qtcf.second.resize(boost::extents[m_block_count][m_block_size][n]);
	}
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate binary correlation functions results");
    }
}

/**
 * returns simulation time belonging to block sample at given block offset
 */
template <unsigned dimension, typename T>
double correlation<dimension, T>::block_sample_time(unsigned int block, unsigned int offset) const
{
    double time = m_timestep * std::pow(m_block_size, block / 2) * (offset + 1);

    if (block % 2) {
	// shifted block
	time *= m_block_shift;
    }
    return time;
}

/**
 * create HDF5 correlations output file
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::open(std::string const& filename)
{
    LOG("write correlations to file: " << filename);
    try {
	// truncate existing file
	m_file = H5::H5File(filename, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 correlations output file");
    }
    // create parameter group
    m_file.createGroup("param");
}

/**
 * close HDF5 file
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::close()
{
    try {
	m_file.close();
    }
    catch (H5::Exception const& e) {
	throw exception("failed to close HDF5 correlations output file");
    }
}

/**
 * returns HDF5 parameter group
 */
template <unsigned dimension, typename T>
H5param correlation<dimension, T>::attrs()
{
    return H5param(m_file.openGroup("param"));
}

/**
 * write parameters to HDF5 parameter group
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::attrs(H5::Group const& param) const
{
    H5xx::group node(param.createGroup("correlation"));
    node["steps"] = m_steps;
    node["time"] = m_time;
    node["timestep"] = m_timestep;
    node["block_size"] = m_block_size;
    node["block_shift"] = m_block_shift;
    node["block_count"] = m_block_count;
    node["max_samples"] = m_max_samples;
    node["q_values"] = m_q_vector.size();
}

/**
 * check if sample is acquired for given simulation step
 */
template <unsigned dimension, typename T>
bool correlation<dimension, T>::sample(uint64_t const& step) const
{
    for (unsigned int i = 0; i < m_block_count; ++i) {
	if (m_block_samples[i] >= m_max_samples)
	    continue;
	if (step % m_block_freq[i])
	    continue;

	return true;
    }
    return false;
}

/**
 * sample time correlation functions
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::sample(std::vector<T> const& r, std::vector<T> const& v, uint64_t const& step)
{
    sample_type sample(r, v, m_q_vector);

    for (unsigned int i = 0; i < m_block_count; ++i) {
	if (m_block_samples[i] >= m_max_samples)
	    continue;
	if (step % m_block_freq[i])
	    continue;

	m_block[i].push_back(sample);

	if (m_block[i].full()) {
	    autocorrelate_block(i);
	    if (i < 2) {
		// sample only full blocks in lowest levels to account for strong correlations
		m_block[i].clear();
	    }
	    m_block_samples[i]++;
	    if (m_max_samples == m_block_samples[i]) {
		LOG("finished sampling on block level " << i << " at step " << step);
		// schedule remaining MD simulation runtime estimate
		alarm(60);
	    }
	}
    }
}

/**
 * apply correlation functions to block samples
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::autocorrelate_block(unsigned int n)
{
    foreach (tcf_pair& tcf, m_tcf) {
	boost::apply_visitor(tcf_apply_visitor_gen(m_block[n].begin(), m_block[n].end(), tcf.second[n].begin()), tcf.first);
    }
    foreach (qtcf_pair& qtcf, m_qtcf) {
	boost::apply_visitor(tcf_apply_visitor_gen(m_block[n].begin(), m_block[n].end(), qtcf.second[n].begin()), qtcf.first);
    }
}

/**
 * write correlation function results to HDF5 file
 */
template <unsigned dimension, typename T>
void correlation<dimension, T>::write()
{
    // compute higher block correlations for remaining samples
    for (unsigned int i = 2; i < m_block_count; ++i) {
	while (m_block[i].size() > 2) {
	    m_block[i].pop_front();
	    autocorrelate_block(i);
	}
    }

    try {
	// assure adequate number of samples per block
	unsigned int max_blocks;
	for (max_blocks = 0; max_blocks < m_block.size(); ++max_blocks) {
	    if (m_block_samples[max_blocks] < 1) {
		LOG_WARNING("could gather only " << max_blocks << " blocks of correlation function results");
		break;
	    }
	}

	if (max_blocks == 0)
	    return;

	H5::DataType tid(H5::PredType::NATIVE_DOUBLE);

	// iterate over correlation functions
	foreach (tcf_pair& tcf, m_tcf) {
	    // dataspace for correlation function results
	    hsize_t dim[3] = { max_blocks, tcf.second.shape()[1], 3 };
	    H5::DataSpace ds(3, dim);
	    // correlation function name
	    char const* name = boost::apply_visitor(tcf_name_visitor(), tcf.first);

	    // create dataset for correlation function results
	    H5::DataSet dataset(m_file.createDataSet(name, tid, ds));

	    // compose results in memory
	    boost::multi_array<double, 3> data(boost::extents[dim[0]][dim[1]][dim[2]]);

	    for (unsigned int j = 0; j < data.size(); ++j) {
		for (unsigned int k = 0; k < data[j].size(); ++k) {
		    // time interval
		    data[j][k][0] = block_sample_time(j, k);
		    // mean average
		    data[j][k][1] = tcf.second[j][k].mean();
		    // standard error of mean
		    data[j][k][2] = tcf.second[j][k].err();
		}
	    }

	    // write results to HDF5 file
	    dataset.write(data.data(), tid);
	}

	// iterate over binary correlation functions
	foreach (qtcf_pair& qtcf, m_qtcf) {
	    // dataspace for binary correlation function results
	    hsize_t dim[4] = { m_q_vector.size(), max_blocks, qtcf.second.shape()[1], 4 };
	    H5::DataSpace ds(4, dim);
	    // correlation function name
	    char const* name = boost::apply_visitor(tcf_name_visitor(), qtcf.first);

	    // create dataset for correlation function results
	    H5::DataSet dataset(m_file.createDataSet(name, tid, ds));

	    // compose results in memory
	    boost::multi_array<double, 4> data(boost::extents[dim[0]][dim[1]][dim[2]][dim[3]]);

	    for (unsigned int j = 0; j < data.size(); ++j) {
		for (unsigned int k = 0; k < data[j].size(); ++k) {
		    for (unsigned int l = 0; l < data[j][k].size(); ++l) {
			// q-value
			data[j][k][l][0] = m_q_vector[j];
			// time interval
			data[j][k][l][1] = (l > 0) ? block_sample_time(k, l - 1) : 0;
			// mean average
			data[j][k][l][2] = qtcf.second[k][l][j].mean();
			// standard error of mean
			data[j][k][l][3] = qtcf.second[k][l][j].err();
		    }
		}
	    }

	    // write results to HDF5 file
	    dataset.write(data.data(), tid);
	}
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write results to correlations file");
    }
}

} // namespace mdsim

#undef foreach

#endif /* ! MDSIM_AUTOCORRELATION_HPP */
