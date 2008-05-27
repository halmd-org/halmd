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
#include "exception.hpp"
#include "log.hpp"
#include "tcf.hpp"


namespace mdsim {

template <typename T>
struct phase_space_point
{
    typedef T vector_type;

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
    typedef phase_space_point<std::vector<T> > phase_space_type;
    typedef phase_space_samples<std::vector<T> > block_type;
    typedef typename std::vector<block_type>::iterator block_iterator;
    typedef typename std::vector<block_type>::const_iterator block_const_iterator;

    typedef typename boost::make_variant_over<tcf_types>::type tcf_type;
    typedef typename boost::array<tcf_type, 3> tcf_array;
    typedef typename tcf_array::iterator tcf_iterator;
    typedef typename tcf_array::const_iterator tcf_const_iterator;

    typedef typename boost::multi_array<accumulator<double>, 2> result_type;
    typedef typename boost::array<result_type, 3> result_array;
    typedef typename result_array::iterator result_iterator;
    typedef typename result_array::const_iterator result_const_iterator;

public:
    /** compute block parameters */
    autocorrelation(H5param const& param);
    /** create HDF5 correlations output file */
    void open(std::string const& filename);
    /** close HDF5 correlations output file */
    void close();

    /** get number of simulation steps */
    uint64_t steps() const;
    /** get block size */
    unsigned int block_size() const;
    /** get block shift */
    unsigned int block_shift() const;
    /** get block count */
    unsigned int block_count() const;
    /** get maximum number of samples per block */
    uint64_t max_samples() const;

    /** write global simulation parameters to autocorrelation output file */
    void write(H5param const& param);

    void sample(std::vector<T> const& r, std::vector<T> const& v);
    void write();

private:
    void autocorrelate(phase_space_type const& sample, unsigned int offset);
    void autocorrelate_block(unsigned int n);
    void finalize();
    void compute_block_param(unsigned int block_size_, uint64_t steps, uint64_t max_samples_);
    double timegrid(unsigned int n, unsigned int k);

private:
    /** phase space sample blocks */
    std::vector<block_type> block;
    /** correlation functions results */
    result_array result;
    /** correlation functions */
    tcf_array tcf;
    /** HDF5 output file */
    H5::H5File file;

    /** number of simulation steps */
    uint64_t steps_;
    /** block size */
    unsigned int block_size_;
    /** block shift */
    unsigned int block_shift_;
    /** block count */
    unsigned int block_count_;
    /** maximum number of samples per block */
    uint64_t max_samples_;
    /** simulation timestep */
    double timestep_;
};


/**
 * compute block parameters
 */
template <unsigned dimension, typename T>
autocorrelation<dimension, T>::autocorrelation(H5param const& param) : steps_(param.steps()), block_size_(param.block_size()), max_samples_(param.max_samples()), timestep_(param.timestep())
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    // compute block shift
    block_shift_ = std::floor(std::sqrt(block_size_));

    // compute block level count
    block_count_ = 0;
    if (block_size_ > 1) {
	for (unsigned int n = block_size_; n <= steps_; n *= block_size_) {
	    ++block_count_;
	    if ((n * block_shift_) > steps_)
		break;
	    ++block_count_;
	}
    }

    LOG("block size  = " << block_size_);
    LOG("block shift = " << block_shift_);
    LOG("block count = " << block_count_);
    LOG("max samples = " << max_samples_);

    // validate block parameters
    if (max_samples_ < block_size_) {
	throw exception("maximum number of samples must not be smaller than block size");
    }
    if (block_shift_ < 2) {
	throw exception("computed block shift is less than 2, larger block size required");
    }
    if (block_count_ < 2) {
	throw exception("computed block count is less than 2, more simulations steps required");
    }

    // allocate phase space sample blocks
    try {
	block.resize(block_count_, block_type(block_size_));
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
	    it->resize(boost::extents[block_count_][block_size_ - 1]);
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
 * get number of simulation steps
 */
template <unsigned dimension, typename T>
uint64_t autocorrelation<dimension, T>::steps() const
{
    return steps_;
}

/**
 * get number of simulation steps
 */
template <unsigned dimension, typename T>
unsigned int autocorrelation<dimension, T>::block_size() const
{
    return block_size_;
}

/**
 * get block shift
 */
template <unsigned dimension, typename T>
unsigned int autocorrelation<dimension, T>::block_shift() const
{
    return block_shift_;
}

/**
 * get block count
 */
template <unsigned dimension, typename T>
unsigned int autocorrelation<dimension, T>::block_count() const
{
    return block_count_;
}

/**
 * get maximum number of samples per block
 */
template <unsigned dimension, typename T>
uint64_t autocorrelation<dimension, T>::max_samples() const
{
    return max_samples_;
}

/**
 * write global simulation parameters to autocorrelation output file
 */
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::write(H5param const& param)
{
    param.write(file.createGroup("/parameters"));
}

template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::sample(std::vector<T> const& r, std::vector<T> const& v)
{
    // sample odd level blocks
    autocorrelate(phase_space_type(r, v), 0);

    if (0 == block[0].count % block_shift_) {
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
    if ((0 == block[offset].count % block_size_) && block[offset].nsample < max_samples_) {
	autocorrelate_block(offset);
	block[offset].nsample++;
    }

    for (unsigned int i = offset + 2; i < block_count_; i += 2) {
	// check if coarse graining is possible
	if (block[i - 2].count % block_size_) {
	    break;
	}

	// add phase space sample from lower level block middle
	block[i].samples.push_back(block[i - 2].samples[block_size_ / 2]);
	block[i].count++;

	// autocorrelate block if circular buffer is full
	if (block[i].samples.full() && block[i].nsample < max_samples_) {
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
    for (unsigned int i = 2; i < block_count_; ++i) {
	while (block[i].nsample < max_samples_ && block[i].samples.size() > 2) {
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


template <unsigned dimension, typename T>
double autocorrelation<dimension, T>::timegrid(unsigned int block, unsigned int sample)
{
    if (block % 2) {
	// shifted block
	return timestep_ * powf(block_size_, block / 2) * (sample + 1) * block_shift_;
    }

    return timestep_ * powf(block_size_, block / 2) * (sample + 1);
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
	    H5::DataSet set(file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds));

	    // compose results in memory
	    boost::multi_array<double, 3> data(boost::extents[dim[0]][dim[1]][dim[2]]);

	    for (unsigned int j = 0; j < data.size(); ++j) {
		for (unsigned int k = 0; k < data[j].size(); ++k) {
		    // time interval
		    data[j][k][0] = timegrid(j, k);
		    // mean average
		    data[j][k][1] = result[i][j][k].mean();
		    // standard error of mean
		    data[j][k][2] = result[i][j][k].err();
		}
	    }

	    // write results to HDF5 file
	    set.write(data.data(), H5::PredType::NATIVE_DOUBLE);
	}

	file.close();
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write results to correlations file");
    }
}

} // namespace mdsim

#endif /* ! MDSIM_AUTOCORRELATION_HPP */
