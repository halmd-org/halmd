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

#include <boost/array.hpp>
// requires patch from http://svn.boost.org/trac/boost/ticket/1852
#include <boost/circular_buffer.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/multi_array.hpp>
#include <boost/variant.hpp>
#include <hdf5.hpp>
#include <string>
#include <vector>
#include "accumulator.hpp"
#include "exception.hpp"
#include "log.hpp"
#include "options.hpp"
#include "tcf.hpp"
#include "trajectory.hpp"


namespace mdsim {

template <typename T>
struct phase_space_samples
{
    /** block samples */
    boost::circular_buffer<phase_space_point<T> > samples;
    /** coarse grain count */
    unsigned int count;
    /** sample count */
    unsigned int nsample;

    phase_space_samples(size_t size) : samples(size), count(0), nsample(0) { }
};


template <unsigned dimension, typename S>
class autocorrelation
{
private:
    typedef phase_space_point<typename S::vector_type> phase_space_type;
    typedef phase_space_samples<typename S::vector_type> block_type;
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
    autocorrelation(options const& opts);

    /** get number of simulation steps */
    uint64_t steps() const;

    /** read autocorrelation parameters from HDF5 file */
    void read_param(H5::Group const& root);
    /** write autocorrelation parameters to HDF5 file */
    void write_param(H5::Group& root) const;

    void sample(S const& s);
    void write(std::string const& path, float timestep);

private:
    void sample(S const& s, unsigned int offset);
    void autocorrelate_block(unsigned int n);
    void finalize();
    void compute_block_param(unsigned int block_size_, uint64_t steps, uint64_t max_samples_);
    float timegrid(unsigned int n, unsigned int k, float timestep);

private:
    /** phase space sample blocks */
    std::vector<block_type> block;
    /** correlation functions results */
    result_array result;
    /** correlation functions */
    tcf_array tcf;

    /** number of simulation steps */
    uint64_t steps_;
    /** block size */
    unsigned int block_size;
    /** block shift */
    unsigned int block_shift;
    /** block count */
    unsigned int block_count;
    /** maximum number of samples per block */
    uint64_t max_samples;
};


template <unsigned dimension, typename S>
autocorrelation<dimension, S>::autocorrelation(options const& opts)
{
    // compute block parameters
    compute_block_param(opts.block_size(), opts.steps(), opts.max_samples());

    // allocate phase space sample blocks
    try {
	block.resize(block_count, block_type(block_size));
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
	    it->resize(boost::extents[block_count][block_size - 1]);
	}
    }
    catch (std::bad_alloc const& e) {
	throw exception("failed to allocate correlation functions results");
    }
}


/**
 * compute block parameters
 */
template <unsigned dimension, typename S>
void autocorrelation<dimension, S>::compute_block_param(unsigned int block_size_, uint64_t steps, uint64_t max_samples_)
{
    // set number of simulation steps
    steps_ = steps;
    // set block size
    block_size = block_size_;
    // compute block shift
    block_shift = std::floor(std::sqrt(block_size));
    // compute block level count
    block_count = 0;
    if (block_size > 1) {
	for (unsigned int n = block_size; n <= steps_; n *= block_size) {
	    ++block_count;
	    if ((n * block_shift) > steps_)
		break;
	    ++block_count;
	}
    }
    // set maximum number of samples per block
    max_samples = max_samples_;

    LOG("block size  = " << block_size);
    LOG("block shift = " << block_shift);
    LOG("block count = " << block_count);
    LOG("max samples = " << max_samples);

    // validate block parameters
    if (max_samples < block_size) {
	throw exception("maximum number of samples must not be smaller than block size");
    }
    if (block_shift < 2) {
	throw exception("computed block shift is less than 2, larger block size required");
    }
    if (block_count < 2) {
	throw exception("computed block count is less than 2, more simulations steps required");
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
 * read autocorrelation parameters from HDF5 file
 */
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::read_param(H5::Group const& root)
{
    try {
	H5ext::Group param(root.openGroup("autocorrelation"));

	// block parameters
	compute_block_param(param["block_size"].as<unsigned int>(), param["steps"].as<uint64_t>(), param["max_samples"].as<uint64_t>());
    }
    catch (H5::Exception const& e) {
	throw exception("failed to read autocorrelation parameters from HDF5 file");
    }
}

/**
 * write autocorrelation parameters to HDF5 file
 */
template <unsigned dimension, typename T>
void autocorrelation<dimension, T>::write_param(H5::Group& root) const
{
    try {
	H5ext::Group param(root.createGroup("autocorrelation"));

	// number of simulation steps
	param["steps"] = steps_;
	// block size
	param["block_size"] = block_size;
	// block shift
	param["block_shift"] = block_shift;
	// block count
	param["block_count"] = block_count;
	// maximum number of samples per block
	param["max_samples"] = max_samples;
    }
    catch (H5::Exception const& e) {
	throw exception("failed to write autocorrelation parameters to HDF5 file");
    }
}

template <unsigned dimension, typename S>
void autocorrelation<dimension, S>::sample(S const& s)
{
    // sample odd level blocks
    sample(s, 0);

    if (block[0].count % block_shift == 0) {
	// sample even level blocks
	sample(s, 1);
    }
}


template <unsigned dimension, typename S>
void autocorrelation<dimension, S>::sample(S const& s, unsigned int offset)
{
    // add phase space sample to lowest block
    block[offset].samples.push_back(phase_space_type(s.R, s.v));
    block[offset].count++;

    // autocorrelate block if circular buffer has been replaced completely
    if (block[offset].count == block_size && block[offset].nsample < max_samples) {
	autocorrelate_block(offset);
	block[offset].nsample++;
    }

    for (unsigned int i = offset + 2; i < block_count; i += 2) {
	// check if coarse graining is possible
	if (block[i - 2].count < block_size) {
	    break;
	}

	// add phase space sample from lower level block middle
	block[i].samples.push_back(block[i - 2].samples[block_size / 2]);
	block[i].count++;
	// reset lower level block coarse grain count
	block[i - 2].count = 0;

	// autocorrelate block if circular buffer is full
	if (block[i].samples.full() && block[i].nsample < max_samples) {
	    autocorrelate_block(i);
	    block[i].nsample++;
	}
    }
}


/**
 * compute correlations for remaining samples in all blocks
 */
template <unsigned dimension, typename S>
void autocorrelation<dimension, S>::finalize()
{
    for (unsigned int i = 2; i < block_count; ++i) {
	while (block[i].nsample < max_samples && block[i].samples.size() > 2) {
	    block[i].samples.pop_front();
	    autocorrelate_block(i);
	    block[i].nsample++;
	}
    }
}


/**
 * apply correlation functions to block samples
 */
template <unsigned dimension, typename S>
void autocorrelation<dimension, S>::autocorrelate_block(unsigned int n)
{
    for (unsigned int i = 0; i < tcf.size(); ++i) {
	boost::apply_visitor(gen_tcf_apply_visitor(block[n].samples.begin(), block[n].samples.end(), result[i][n].begin()), tcf[i]);
    }
}


template <unsigned dimension, typename S>
float autocorrelation<dimension, S>::timegrid(unsigned int block, unsigned int sample, float timestep)
{
    if (block % 2) {
	// shifted block
	return timestep * powf(block_size, block / 2) * (sample + 1) * block_shift;
    }

    return timestep * powf(block_size, block / 2) * (sample + 1);
}


/**
 * write correlation function results to HDF5 file
 */
template <unsigned dimension, typename S>
void autocorrelation<dimension, S>::write(std::string const& path, float timestep)
{
    // compute correlations for remaining samples in all blocks
    finalize();

#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    H5::H5File file;

    try {
	// create empty HDF5 file
	file = H5::H5File(path, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create correlations file");
    }

    try {
	H5ext::Group root(file.openGroup("/"));

	// write block parameters
	root["block_count"] = block_count;
	root["block_size"] = block_size;
	root["block_shift"] = block_shift;
	root["max_samples"] = max_samples;

	// write simulation parameters
	root["timestep"] = timestep;
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write attributes to correlations file");
    }

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
		    data[j][k][0] = timegrid(j, k, timestep);
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
