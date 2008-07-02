/* Time correlation functions
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

#ifndef MDSIM_AUTOCORRELATION_HPP
#define MDSIM_AUTOCORRELATION_HPP

#include <H5Cpp.h>
#include <boost/array.hpp>
// requires patch from http://svn.boost.org/trac/boost/ticket/1852
#include <boost/circular_buffer.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/variant.hpp>
#include <cmath>
#include <string>
#include <vector>
#include "H5param.hpp"
#include "H5xx.hpp"
#include "accumulator.hpp"
#include "block.hpp"
#include "exception.hpp"
#include "log.hpp"
#include "tcf.hpp"


#define foreach BOOST_FOREACH

namespace mdsim {

/**
 * Phase space sample
 */
template <unsigned dimension, typename T, typename U>
struct phase_space_point
{
    // swappable host memory vector type
    typedef std::vector<T> vector_type;
    typedef typename T::value_type value_type;
    typedef typename std::vector<std::pair<T, T> > density_vector_type;
    typedef typename std::vector<value_type> q_values_type;
    typedef typename std::vector<value_type>::const_iterator q_values_const_iterator;
    typedef typename std::pair<q_values_const_iterator, q_values_const_iterator> q_values_const_iterator_pair;

    phase_space_point(cuda::host::vector<U> const& h_r, cuda::host::vector<U> const& h_v, q_values_type const& q) : q(q.begin(), q.end()), rho(q.size(), std::pair<T, T>(0, 0))
    {
	r.reserve(h_r.size());
	v.reserve(h_v.size());

	for (size_t i = 0; i < h_r.size(); ++i) {
	    // convert from GPU type to host type
	    const T r(h_r[i]), v(h_v[i]);

	    this->r.push_back(r);
	    this->v.push_back(v);

	    // spatial Fourier transformation
	    for (unsigned int j = 0; j < q.size(); ++j) {
		// compute averages to maintain accuracy with single precision floating-point
		rho[j].first += (cos(r * q[j]) - rho[j].first) / (i + 1);
		rho[j].second += (sin(r * q[j]) - rho[j].second) / (i + 1);
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
    /** q-values */
    q_values_const_iterator_pair q;
    /** spatially Fourier transformed density for given q-values */
    density_vector_type rho;
};

/**
 * Block of phase space samples
 */
template <unsigned dimension, typename T, typename U>
struct phase_space_samples : boost::circular_buffer<phase_space_point<dimension, T, U> >
{
    phase_space_samples(size_t size) : boost::circular_buffer<phase_space_point<dimension, T, U> >(size), count(0), samples(0) { }

    /** trajectory sample count */
    uint64_t count;
    /** block autocorrelation count */
    uint64_t samples;
};

/**
 * Autocorrelation block algorithm
 */
template <unsigned dimension, typename T, typename U>
class autocorrelation
{
public:
    /** sample vector in page-locked host memory */
    typedef cuda::host::vector<U> vector_type;

    /** phase space sample type */
    typedef phase_space_point<dimension, T, U> phase_space_type;
    /** phase space samples block type */
    typedef phase_space_samples<dimension, T, U> block_type;

    /** correlation function type */
    typedef typename boost::make_variant_over<tcf_types>::type tcf_type;
    /** correlation function results type */
    typedef boost::multi_array<accumulator<float>, 2> tcf_result_type;
    /** correlation function and results pair type */
    typedef std::pair<tcf_type, tcf_result_type> tcf_pair;

    /** binary correlation function type */
    typedef typename boost::make_variant_over<qtcf_types>::type qtcf_type;
    /** binary correlation function results type */
    typedef boost::multi_array<accumulator<float>, 3> qtcf_result_type;
    /** binary correlation function and results pair type */
    typedef std::pair<qtcf_type, qtcf_result_type> qtcf_pair;

public:
    /** initialize correlation functions */
    autocorrelation(block_param<dimension, T> const& param, float const& box, unsigned int nq);

    /** create HDF5 correlations output file */
    void open(std::string const& filename);
    /** returns HDF5 parameter group */
    H5param attrs();
    /** write parameters to HDF5 parameter group */
    void attrs(H5::Group const& param) const;
    /** sample time correlation functions */
    void sample(vector_type const& r, vector_type const& v);
    /** write correlation function results to HDF5 file */
    void write();
    /** close HDF5 file */
    void close();


private:
    /** autocorrelate odd or even blocks */
    void autocorrelate(phase_space_type const& sample, unsigned int offset);
    /** apply correlation functions to block samples */
    void autocorrelate_block(unsigned int n);
    /** compute correlations for remaining samples in all blocks */
    void finalize();

private:
    /** block algorithm parameters */
    block_param<dimension, T> param;
    /** q-values for spatial Fourier transformation */
    std::vector<float> q_;

    /** phase space sample blocks */
    std::vector<block_type> block;
    /** correlation functions and results */
    boost::array<tcf_pair, 3> tcf_;
    /** binary correlation functions and results */
    boost::array<qtcf_pair, 2> qtcf_;

    /** HDF5 output file */
    H5::H5File file;
};

/**
 * initialize correlation functions
 */
template <unsigned dimension, typename T, typename U>
autocorrelation<dimension, T, U>::autocorrelation(block_param<dimension, T> const& param, float const& box, unsigned int nq) : param(param)
{
    // allocate phase space sample blocks
    try {
	block.resize(param.block_count(), block_type(param.block_size()));
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate phase space sample blocks");
    }

    // compute q-values for spatial Fourier transformation */
    for (unsigned int k = 1; k <= nq; ++k) {
	// integer multiple of smallest q-value
	q_.push_back(k * 2 * M_PI / box);
    }

    // setup correlation functions
    tcf_[0].first = mean_square_displacement();
    tcf_[1].first = mean_quartic_displacement();
    tcf_[2].first = velocity_autocorrelation();
    // setup binary correlation functions
    qtcf_[0].first = intermediate_scattering_function();
    qtcf_[1].first = self_intermediate_scattering_function();

    try {
	// allocate correlation functions results
	foreach (tcf_pair& tcf, tcf_) {
	    tcf.second.resize(boost::extents[param.block_count()][param.block_size() - 1]);
	}
	foreach (qtcf_pair& qtcf, qtcf_) {
	    qtcf.second.resize(boost::extents[param.block_count()][param.block_size()][q_.size()]);
	}
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate correlation functions results");
    }
}

/**
 * create HDF5 correlations output file
 */
template <unsigned dimension, typename T, typename U>
void autocorrelation<dimension, T, U>::open(std::string const& filename)
{
    LOG("write correlations to file: " << filename);
    try {
	// truncate existing file
	file = H5::H5File(filename, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 correlations output file");
    }
    // create parameter group
    file.createGroup("param");
}

/**
 * returns HDF5 parameter group
 */
template <unsigned dimension, typename T, typename U>
H5param autocorrelation<dimension, T, U>::attrs()
{
    return H5param(file.openGroup("param"));
}

/**
 * write parameters to HDF5 parameter group
 */
template <unsigned dimension, typename T, typename U>
void autocorrelation<dimension, T, U>::attrs(H5::Group const& param) const
{
    H5xx::group node(param.openGroup("autocorrelation"));
    node["q_values"] = q_.size();
}

/**
 * sample time correlation functions
 */
template <unsigned dimension, typename T, typename U>
void autocorrelation<dimension, T, U>::sample(vector_type const& r, vector_type const& v)
{
    // sample odd level blocks
    autocorrelate(phase_space_type(r, v, q_), 0);

    if (0 == block[0].count % param.block_shift()) {
	// sample even level blocks
	autocorrelate(phase_space_type(r, v, q_), 1);
    }
}


/**
 * autocorrelate odd or even blocks
 */
template <unsigned dimension, typename T, typename U>
void autocorrelation<dimension, T, U>::autocorrelate(phase_space_type const& sample, unsigned int offset)
{
    // add phase space sample to lowest block
    block[offset].push_back(sample);
    block[offset].count++;

    // autocorrelate block if circular buffer has been replaced completely
    if ((0 == block[offset].count % param.block_size()) && block[offset].samples < param.max_samples()) {
	autocorrelate_block(offset);
	block[offset].samples++;
    }

    for (unsigned int i = offset + 2; i < param.block_count(); i += 2) {
	// check if coarse graining is possible
	if (block[i - 2].count % param.block_size()) {
	    break;
	}

	// add phase space sample from lower level block middle
	block[i].push_back(block[i - 2][param.block_size() / 2]);
	block[i].count++;

	// autocorrelate block if circular buffer is full
	if (block[i].full() && block[i].samples < param.max_samples()) {
	    autocorrelate_block(i);
	    block[i].samples++;
	}
    }
}

/**
 * compute correlations for remaining samples in all blocks
 */
template <unsigned dimension, typename T, typename U>
void autocorrelation<dimension, T, U>::finalize()
{
    for (unsigned int i = 2; i < param.block_count(); ++i) {
	while (block[i].samples < param.max_samples() && block[i].size() > 2) {
	    block[i].pop_front();
	    autocorrelate_block(i);
	}
    }
}

/**
 * apply correlation functions to block samples
 */
template <unsigned dimension, typename T, typename U>
void autocorrelation<dimension, T, U>::autocorrelate_block(unsigned int n)
{
    foreach (tcf_pair& tcf, tcf_) {
	boost::apply_visitor(tcf_apply_visitor_gen(block[n].begin(), block[n].end(), tcf.second[n].begin()), tcf.first);
    }
    foreach (qtcf_pair& qtcf, qtcf_) {
	boost::apply_visitor(tcf_apply_visitor_gen(block[n].begin(), block[n].end(), qtcf.second[n].begin()), qtcf.first);
    }
}

/**
 * write correlation function results to HDF5 file
 */
template <unsigned dimension, typename T, typename U>
void autocorrelation<dimension, T, U>::write()
{
    // compute correlations for remaining samples in all blocks
    finalize();

    try {
	// assure adequate number of samples per block
	unsigned int max_blocks;
	for (max_blocks = 0; max_blocks < block.size(); ++max_blocks) {
	    if (block[max_blocks].samples < 1) {
		LOG_WARNING("could gather only " << max_blocks << " blocks of correlation function results");
		break;
	    }
	}

	if (max_blocks == 0)
	    return;

	// iterate over correlation functions
	foreach (tcf_pair& tcf, tcf_) {
	    // dataspace for correlation function results
	    hsize_t dim[3] = { max_blocks, tcf.second.shape()[1], 3 };
	    H5::DataSpace ds(3, dim);
	    // correlation function name
	    char const* name = boost::apply_visitor(tcf_name_visitor(), tcf.first);

	    // create dataset for correlation function results
	    H5::DataSet dataset(file.createDataSet(name, H5::PredType::NATIVE_FLOAT, ds));

	    // compose results in memory
	    boost::multi_array<float, 3> data(boost::extents[dim[0]][dim[1]][dim[2]]);

	    for (unsigned int j = 0; j < data.size(); ++j) {
		for (unsigned int k = 0; k < data[j].size(); ++k) {
		    // time interval
		    data[j][k][0] = param.timegrid(j, k);
		    // mean average
		    data[j][k][1] = tcf.second[j][k].mean();
		    // standard error of mean
		    data[j][k][2] = tcf.second[j][k].err();
		}
	    }

	    // write results to HDF5 file
	    dataset.write(data.data(), H5::PredType::NATIVE_FLOAT);
	}

	// iterate over binary correlation functions
	foreach (qtcf_pair& qtcf, qtcf_) {
	    // dataspace for binary correlation function results
	    hsize_t dim[4] = { q_.size(), max_blocks, qtcf.second.shape()[1], 4 };
	    H5::DataSpace ds(4, dim);
	    // correlation function name
	    char const* name = boost::apply_visitor(tcf_name_visitor(), qtcf.first);

	    // create dataset for correlation function results
	    H5::DataSet dataset(file.createDataSet(name, H5::PredType::NATIVE_FLOAT, ds));

	    // compose results in memory
	    boost::multi_array<float, 4> data(boost::extents[dim[0]][dim[1]][dim[2]][dim[3]]);

	    for (unsigned int j = 0; j < data.size(); ++j) {
		for (unsigned int k = 0; k < data[j].size(); ++k) {
		    for (unsigned int l = 0; l < data[j][k].size(); ++l) {
			// q-value
			data[j][k][l][0] = q_[j];
			// time interval
			data[j][k][l][1] = (l > 0) ? param.timegrid(k, l - 1) : 0;
			// mean average
			data[j][k][l][2] = qtcf.second[k][l][j].mean();
			// standard error of mean
			data[j][k][l][3] = qtcf.second[k][l][j].err();
		    }
		}
	    }

	    // write results to HDF5 file
	    dataset.write(data.data(), H5::PredType::NATIVE_FLOAT);
	}
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write results to correlations file");
    }
}

/**
 * close HDF5 file
 */
template <unsigned dimension, typename T, typename U>
void autocorrelation<dimension, T, U>::close()
{
    try {
	file.close();
    }
    catch (H5::Exception const& e) {
	throw exception("failed to close HDF5 correlations output file");
    }
}

} // namespace mdsim

#undef foreach

#endif /* ! MDSIM_AUTOCORRELATION_HPP */
