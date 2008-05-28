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
#include <boost/mpl/vector.hpp>
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


namespace mdsim {

template <typename T>
struct phase_space_point
{
    typedef T vector_type;

    phase_space_point() {}
    phase_space_point(T const& r, T const& v) : r(r), v(v) {}

    /** particle positions */
    T r;
    /** particle velocities */
    T v;
};


template <typename T>
struct phase_space_samples
{
    phase_space_samples(size_t size) : samples(size), count(0), nsample(0) { }

    /** block samples */
    boost::circular_buffer<phase_space_point<T> > samples;
    /** trajectory sample count */
    unsigned int count;
    /** block autocorrelation count */
    unsigned int nsample;
};


template <unsigned dimension, typename T>
class autocorrelation
{
private:
    typedef phase_space_point<cuda::host::vector<T> > phase_space_type;
    typedef phase_space_samples<cuda::host::vector<T> > block_type;
    typedef typename std::vector<block_type>::iterator block_iterator;
    typedef typename std::vector<block_type>::const_iterator block_const_iterator;

    typedef typename boost::make_variant_over<tcf_types>::type tcf_type;
    typedef typename boost::array<tcf_type, 3> tcf_array;
    typedef typename tcf_array::iterator tcf_iterator;
    typedef typename tcf_array::const_iterator tcf_const_iterator;

    typedef typename boost::multi_array<accumulator<float>, 2> result_type;
    typedef typename boost::array<result_type, 3> result_array;
    typedef typename result_array::iterator result_iterator;
    typedef typename result_array::const_iterator result_const_iterator;

public:
    /** initialize correlation functions */
    autocorrelation(block_param<dimension, T> const& param);

    /** create HDF5 correlations output file */
    void open(std::string const& filename);
    /** close HDF5 correlations output file */
    void close();

    /** dump global simulation parameters to HDF5 file */
    autocorrelation<dimension, T>& operator<<(H5param const& param);

    void sample(cuda::host::vector<T> const& r, cuda::host::vector<T> const& v);
    void write();

private:
    void autocorrelate(phase_space_type const& sample, unsigned int offset);
    void autocorrelate_block(unsigned int n);
    void finalize();
    float timegrid(unsigned int n, unsigned int k);

private:
    /** phase space sample blocks */
    std::vector<block_type> block;
    /** correlation functions results */
    result_array result;
    /** correlation functions */
    tcf_array tcf;
    /** HDF5 output file */
    H5::H5File file;

    /** block algorithm parameters */
    block_param<dimension, T> param;
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
    tcf[0] = mean_square_displacement();
    tcf[1] = mean_quartic_displacement();
    tcf[2] = velocity_autocorrelation();

    // allocate correlation functions results
    try {
	for (result_iterator it = result.begin(); it != result.end(); ++it) {
	    it->resize(boost::extents[param.block_count()][param.block_size() - 1]);
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
 * close HDF5 correlations output file
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

/**
 * dump global simulation parameters to HDF5 file
 */
template <unsigned dimension, typename T>
autocorrelation<dimension, T>& autocorrelation<dimension, T>::operator<<(H5param const& param)
{
    param.write(file.createGroup("/parameters"));
    return *this;
}

template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::sample(cuda::host::vector<T> const& r, cuda::host::vector<T> const& v)
{
    // sample odd level blocks
    autocorrelate(phase_space_type(r, v), 0);

    if (0 == block[0].count % param.block_shift()) {
	// sample even level blocks
	autocorrelate(phase_space_type(r, v), 1);
    }
}


template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::autocorrelate(phase_space_type const& sample, unsigned int offset)
{
    // add phase space sample to lowest block
    block[offset].samples.push_back(sample);
    block[offset].count++;

    // autocorrelate block if circular buffer has been replaced completely
    if ((0 == block[offset].count % param.block_size()) && block[offset].nsample < param.max_samples()) {
	autocorrelate_block(offset);
	block[offset].nsample++;
    }

    for (unsigned int i = offset + 2; i < param.block_count(); i += 2) {
	// check if coarse graining is possible
	if (block[i - 2].count % param.block_size()) {
	    break;
	}

	// add phase space sample from lower level block middle
	block[i].samples.push_back(block[i - 2].samples[param.block_size() / 2]);
	block[i].count++;

	// autocorrelate block if circular buffer is full
	if (block[i].samples.full() && block[i].nsample < param.max_samples()) {
	    autocorrelate_block(i);
	    block[i].nsample++;
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
	while (block[i].nsample < param.max_samples() && block[i].samples.size() > 2) {
	    block[i].samples.pop_front();
	    autocorrelate_block(i);
	    block[i].nsample++;
	}
    }
}


/**
 * apply correlation functions to block samples
 */
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::autocorrelate_block(unsigned int n)
{
    for (unsigned int i = 0; i < tcf.size(); ++i) {
	boost::apply_visitor(gen_tcf_apply_visitor(block[n].samples.begin(), block[n].samples.end(), result[i][n].begin()), tcf[i]);
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
	// iterate over correlation functions
	for (unsigned int i = 0; i < tcf.size(); ++i) {
	    // dataspace for correlation function results
	    hsize_t dim[3] = { result[i].shape()[0], result[i].shape()[1], 3 };
	    H5::DataSpace ds(3, dim);
	    // correlation function name
	    char const* name = boost::apply_visitor(tcf_name_visitor(), tcf[i]);

	    // create dataset for correlation function results
	    H5::DataSet set(file.createDataSet(name, H5::PredType::NATIVE_FLOAT, ds));

	    // compose results in memory
	    boost::multi_array<float, 3> data(boost::extents[dim[0]][dim[1]][dim[2]]);

	    for (unsigned int j = 0; j < data.size(); ++j) {
		for (unsigned int k = 0; k < data[j].size(); ++k) {
		    // time interval
		    data[j][k][0] = param.timegrid(j, k);
		    // mean average
		    data[j][k][1] = result[i][j][k].mean();
		    // standard error of mean
		    data[j][k][2] = result[i][j][k].err();
		}
	    }

	    // write results to HDF5 file
	    set.write(data.data(), H5::PredType::NATIVE_FLOAT);
	}

	file.close();
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write results to correlations file");
    }
}

} // namespace mdsim

#endif /* ! MDSIM_AUTOCORRELATION_HPP */
