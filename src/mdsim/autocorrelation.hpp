/* Autocorrelation block algorithm
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
#include <string>
#include <vector>
#include "H5param.hpp"
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
template <typename T>
struct phase_space_point
{
    typedef std::vector<T> vector_type;

    phase_space_point(vector_type const& r, vector_type const& v) : r(r), v(v) {}

    /** particle positions */
    vector_type r;
    /** particle velocities */
    vector_type v;
};

/**
 * Block of phase space samples
 */
template <typename T>
struct phase_space_samples : boost::circular_buffer<phase_space_point<T> >
{
    phase_space_samples(size_t size) : boost::circular_buffer<phase_space_point<T> >(size), count(0), samples(0) { }

    /** trajectory sample count */
    uint64_t count;
    /** block autocorrelation count */
    uint64_t samples;
};

/**
 * Autocorrelation block algorithm
 */
template <unsigned dimension, typename T>
class autocorrelation
{
public:
    /** phase space sample type */
    typedef phase_space_point<T> phase_space_type;
    /** phase space samples block type */
    typedef phase_space_samples<T> block_type;

    /** generic correlation function type */
    typedef typename boost::make_variant_over<tcf_types>::type tcf_type;
    /** correlation function results type */
    typedef boost::multi_array<accumulator<double>, 2> tcf_result_type;
    /** correlation function and results pair type */
    typedef std::pair<tcf_type, tcf_result_type> tcf_pair;

public:
    /** initialize correlation functions */
    autocorrelation(block_param<dimension, T> const& param);

    /** create HDF5 correlations output file */
    void open(std::string const& filename);
    /** dump global simulation parameters to HDF5 file */
    autocorrelation<dimension, T>& operator<<(H5param const& param);
    /** sample trajectory correlation functions */
    void sample(std::vector<T> const& r, std::vector<T> const& v);
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

    /** phase space sample blocks */
    std::vector<block_type> block;
    /** correlation functions and results */
    boost::array<tcf_pair, 3> tcf_;

    /** HDF5 output file */
    H5::H5File file;

};

/**
 * initialize correlation functions
 */
template <unsigned dimension, typename T>
autocorrelation<dimension, T>::autocorrelation(block_param<dimension, T> const& param) : param(param)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    // allocate phase space sample blocks
    try {
	block.resize(param.block_count(), block_type(param.block_size()));
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate phase space sample blocks");
    }

    // setup correlation functions
    tcf_[0].first = mean_square_displacement();
    tcf_[1].first = mean_quartic_displacement();
    tcf_[2].first = velocity_autocorrelation();

    // allocate correlation functions results
    try {
	foreach (tcf_pair& tcf, tcf_) {
	    tcf.second.resize(boost::extents[param.block_count()][param.block_size() - 1]);
	}
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate correlation functions results");
    }
}

/**
 * create HDF5 correlations output file
 */
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::open(std::string const& filename)
{
    LOG("write correlations to file: " << filename);
    try {
	// truncate existing file
	file = H5::H5File(filename, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create HDF5 correlations output file");
    }
}

/**
 * dump global simulation parameters to HDF5 file
 */
template <unsigned dimension, typename T>
autocorrelation<dimension, T>& autocorrelation<dimension, T>::operator<<(H5param const& param)
{
    param.write(file.createGroup("/parameters"));
    return *this;
}

/**
 * sample trajectory correlation functions
 */
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::sample(std::vector<T> const& r, std::vector<T> const& v)
{
    // sample odd level blocks
    autocorrelate(phase_space_type(r, v), 0);

    if (0 == block[0].count % param.block_shift()) {
	// sample even level blocks
	autocorrelate(phase_space_type(r, v), 1);
    }
}


/**
 * autocorrelate odd or even blocks
 */
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::autocorrelate(phase_space_type const& sample, unsigned int offset)
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
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::finalize()
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
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::autocorrelate_block(unsigned int n)
{
    foreach (tcf_pair& tcf, tcf_) {
	boost::apply_visitor(tcf_apply_visitor_gen(block[n].begin(), block[n].end(), tcf.second[n].begin()), tcf.first);
    }
}

/**
 * write correlation function results to HDF5 file
 */
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::write()
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
	    H5::DataSet dataset(file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds));

	    // compose results in memory
	    boost::multi_array<double, 3> data(boost::extents[dim[0]][dim[1]][dim[2]]);

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
	    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
	}
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write results to correlations file");
    }
}

/**
 * close HDF5 file
 */
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::close()
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
