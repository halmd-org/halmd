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

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/for_each.hpp>
#include <ljgpu/sample/correlation.hpp>
#include <ljgpu/sample/tcf_visitor.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <ljgpu/util/exception.hpp>

#define foreach BOOST_FOREACH

namespace ljgpu {

/**
 * set total number of simulation steps
 */
template <int dimension>
void correlation<dimension>::steps(uint64_t value, double timestep)
{
    // set total number of simulation steps
    m_steps = value;
    LOG("total number of simulation steps: " << m_steps);
    // set simulation timestep
    m_timestep = timestep;
    // derive total simulation time
    m_time = value * m_timestep;
    LOG("total simulation time: " << m_time);
}

/**
 * set total simulation time
 */
template <int dimension>
void correlation<dimension>::time(double value, double timestep)
{
    // set total simulation time
    m_time = value;
    LOG("total simulation time: " << m_time);
    // set simulation timestep
    m_timestep = timestep;
    // derive total number of simulation steps
    m_steps = round(m_time / m_timestep);
    LOG("total number of simulation steps: " << m_steps);
}

/**
 * set sample rate for lowest block level
 */
template <int dimension>
void correlation<dimension>::sample_rate(unsigned int sample_rate)
{
    m_sample_rate = sample_rate;
    LOG("sample rate for lowest block level: " << m_sample_rate);
}

/**
 * set maximum number of samples per block
 */
template <int dimension>
void correlation<dimension>::max_samples(uint64_t value)
{
    m_max_samples = value;
    LOG("maximum number of samples per block: " << m_max_samples);
}

/**
 * set minimum number of trajectory samples
 */
template <int dimension>
void correlation<dimension>::min_samples(uint64_t value)
{
    m_min_samples = value;
    LOG("minimum number of trajectory samples: " << m_min_samples);
}

/**
 * set block size
 */
template <int dimension>
void correlation<dimension>::block_size(unsigned int value)
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

    // derive block count and sample frequencies from block size and block shift
    m_block_count = 0;
    m_trajectory_block = 0;
    for (uint64_t n = m_sample_rate, m = n * m_block_shift; ; ++m_block_count) {
	if (m_block_count % 2) {
	    if ((m_steps / m) >= m_min_samples) {
		++m_trajectory_block;
	    }
	    if (m > m_steps) {
		break;
	    }
	    m_block_freq.push_back(m);
	    m *= m_block_size;
	}
	else {
	    if ((m_steps / n) >= m_min_samples) {
		++m_trajectory_block;
	    }
	    if (n > m_steps) {
		break;
	    }
	    m_block_freq.push_back(n);
	    n *= m_block_size;
	}
    }
    LOG("block count: " << m_block_count);
    if (!m_block_count) {
	throw exception("computed block count is zero, more simulations steps required");
    }
    if (!m_trajectory_block) {
	LOG_WARNING("minimum number of trajectory samples not satisfied");
    }
    // convert trajectory block count to block level
    m_trajectory_block = std::max(m_trajectory_block, 1U) - 1;

    m_block_samples.resize(m_block_count, 0);

    // compute block time intervals
    m_block_time.resize(boost::extents[m_block_count][m_block_size]);
    for (unsigned int i = 0; i < m_block_count; ++i) {
	for (unsigned int j = 0; j < m_block_size; ++j) {
	    m_block_time[i][j] = m_timestep * m_sample_rate * std::pow(m_block_size, i / 2) * j;
	    // shifted block
	    if (i % 2) {
		m_block_time[i][j] *= m_block_shift;
	    }
	}
    }
}

/**
 * set q-vectors for spatial Fourier transformation
 */
template <int dimension>
void correlation<dimension>::q_values(std::vector<float> const& values, float error, float box)
{
    typedef std::pair<int, int> q_pair;
    typedef typename q_vector_vector::value_type q_vector_value;

    m_q_value.resize(values.size());
    m_q_vector.resize(values.size());
    m_q_error = error;

    double const q_lattice = 2 * M_PI / box;
    std::vector<q_pair> qq;
    int q_max = 0;

    for (size_t i = 0; i < values.size(); ++i) {
	// adjust q-value to reciprocal lattice
	m_q_value[i] = round(values[i] / q_lattice);
	// upper and lower boundaries within given error
	int q_upper = std::floor(m_q_value[i] * (1 + m_q_error));
	int q_lower = std::ceil(m_q_value[i] * (1 - m_q_error));

	qq.push_back(q_pair(q_lower * q_lower, q_upper * q_upper));
	q_max = std::max(q_max, q_upper);
    }

    find_q_vectors(qq, q_max, m_q_vector);

    for (size_t i = 0; i < values.size(); ++i) {
	foreach(vector_type& q, m_q_vector[i]) {
	    q *= q_lattice;
	}
	m_q_value[i] *= q_lattice;

	LOG("|q| = " << m_q_value[i] << " with " << m_q_vector[i].size() << " vectors");
    }
}

#if WITH_CUDA
template <int dimension>
void correlation<dimension>::add_gpu_correlation_functions(size_t types)
{
    LOG("computing correlation functions on GPU");
    for (size_t type = 0; type < types; ++type) {
	tcf_vector tcf_;
	boost::mpl::for_each<tcf_gpu_types>(boost::bind(&tcf_vector::push_back, boost::ref(tcf_), _1));
	foreach (tcf_variant& tcf, tcf_) {
	    boost::apply_visitor(tcf_set_type(type), tcf);
	}
	m_tcf.insert(m_tcf.end(), tcf_.begin(), tcf_.end());
    }
    m_block.resize(m_block_count, gpu_block_type(m_block_size));
    m_sample = gpu_sample_type();
}
#endif /* WITH_CUDA */

template <int dimension>
void correlation<dimension>::add_host_correlation_functions(size_t types)
{
    LOG("computing correlation functions on host");
    for (size_t type = 0; type < types; ++type) {
	tcf_vector tcf_;
	boost::mpl::for_each<tcf_host_types>(boost::bind(&tcf_vector::push_back, boost::ref(tcf_), _1));
	foreach (tcf_variant& tcf, tcf_) {
	    boost::apply_visitor(tcf_set_type(type), tcf);
	}
	m_tcf.insert(m_tcf.end(), tcf_.begin(), tcf_.end());
    }
    m_block.resize(m_block_count, host_block_type(m_block_size));
    m_sample = host_sample_type();
}

/**
 * create HDF5 correlations output file
 */
template <int dimension>
void correlation<dimension>::open(std::string const& filename, size_t types)
{
    LOG("write correlations to file: " << filename);
    try {
	// truncate existing file
	m_file = H5::H5File(filename, H5F_ACC_TRUNC);
	m_is_open = true;
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 correlations output file");
    }

    try {
	foreach (tcf_variant& tcf, m_tcf) {
	    boost::apply_visitor(tcf_allocate_results(m_block_count, m_block_size, m_q_value.size()), tcf);
	}
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate correlation functions results");
    }

    try {
	foreach (tcf_variant& tcf, m_tcf) {
	    boost::apply_visitor(tcf_create_dataset(m_file, types), tcf);
	}
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create correlation function datasets");
    }
}

/**
 * close HDF5 file
 */
template <int dimension>
void correlation<dimension>::close()
{
    // compute higher block correlations for remaining samples
    for (unsigned int i = 2; i < m_block_count; ++i) {
	while (boost::apply_visitor(tcf_block_size(), m_block[i]) > 2) {
	    boost::apply_visitor(tcf_block_pop_front(), m_block[i]);
	    autocorrelate_block(i);
	}
    }

    // write correlation function results to HDF5 file
    flush();

    try {
	m_file.close();
	m_is_open = false;
    }
    catch (H5::Exception const& e) {
	throw exception("failed to close HDF5 correlations output file");
    }
}

/**
 * write parameters to HDF5 parameter group
 */
template <int dimension>
void correlation<dimension>::param(H5::Group const& param) const
{
    H5xx::group node(param.createGroup("correlation"));
    node["steps"] = m_steps;
    node["time"] = m_time;
    node["sample_rate"] = m_sample_rate;
    node["block_size"] = m_block_size;
    node["block_shift"] = m_block_shift;
    node["block_count"] = m_block_count;
    node["max_samples"] = m_max_samples;
    node["min_samples"] = m_min_samples;
    if (!m_q_value.empty()) {
	node["q_values"] = m_q_value;
	node["q_error"] = m_q_error;
    }
}

/**
 * compute lattice points in first octant on surface of 3-dimensional spheres
 */
template <int dimension>
template <typename T>
typename boost::enable_if<boost::is_same<vector<double, 3>, T>, void>::type
correlation<dimension>::find_q_vectors(std::vector<std::pair<int, int> > const& qq, int q_max, std::vector<std::vector<T> >& q)
{
    // FIXME fast algorithm for lattice points on surface of Ewald's sphere
    for (int x = 0; x <= q_max; ++x) {
	int xx = x * x;
	for (int y = 0; y <= q_max; ++y) {
	    int yy = xx + y * y;
	    for (int z = 0; z <= q_max; ++z) {
		int zz = yy + z * z;
		for (size_t i = 0; i < qq.size(); ++i) {
		    if (zz >= qq[i].first && zz <= qq[i].second) {
			q[i].push_back(vector_type(x, y, z));
		    }
		}
	    }
	}
    }
}

/**
 * compute lattice points in first quadrant on surface of 2-dimensional spheres
 */
template <int dimension>
template <typename T>
typename boost::enable_if<boost::is_same<vector<double, 2>, T>, void>::type
correlation<dimension>::find_q_vectors(std::vector<std::pair<int, int> > const& qq, int q_max, std::vector<std::vector<T> >& q)
{
    // FIXME fast algorithm for lattice points on surface of Ewald's sphere
    for (int x = 0; x <= q_max; ++x) {
	int xx = x * x;
	for (int y = 0; y <= q_max; ++y) {
	    int yy = xx + y * y;
	    for (size_t i = 0; i < qq.size(); ++i) {
		if (yy >= qq[i].first && yy <= qq[i].second) {
		    q[i].push_back(vector_type(x, y));
		}
	    }
	}
    }
}

/**
 * apply correlation functions to block samples
 */
template <int dimension>
void correlation<dimension>::autocorrelate_block(unsigned int n)
{
    foreach (tcf_variant& tcf, m_tcf) {
	boost::apply_visitor(tcf_correlate_block(n, m_q_vector), tcf, m_block[n]);
    }
}

/**
 * write correlation function results to HDF5 file
 */
template <int dimension>
void correlation<dimension>::flush()
{
    // find highest block with adequate number of samples
    unsigned int max_blocks = 0;
    for (; max_blocks < m_block_count; ++max_blocks) {
	if (m_block_samples[max_blocks] < 2)
	    break;
    }
    if (max_blocks < 1)
	return;

    try {
	foreach (tcf_variant& tcf, m_tcf) {
	    boost::apply_visitor(tcf_write_results(m_block_time, m_q_value, max_blocks), tcf);
	}
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write correlation function results");
    }

    try {
	// flush file contents to disk
	m_file.flush(H5F_SCOPE_GLOBAL);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to flush HDF5 correlations file");
    }
}

// explicit template instantiation
template class correlation<2>;
template class correlation<3>;

} // namespace ljgpu
