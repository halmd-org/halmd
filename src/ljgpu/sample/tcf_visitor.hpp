/* Time correlation function visitors
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

#ifndef LJGPU_SAMPLE_TCF_VISITOR_HPP
#define LJGPU_SAMPLE_TCF_VISITOR_HPP

#include <algorithm>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/variant.hpp>
#include <ljgpu/sample/tcf.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <vector>

namespace ljgpu {

/**
 * apply correlation function to block of phase space samples
 */
template <typename T1, typename T2>
class tcf_correlate_block : public boost::static_visitor<>
{
public:
    tcf_correlate_block(unsigned int block, T1 const& sample, T2 const& q_vector)
	: block(block), sample(sample), q_vector(q_vector) {}

    template <typename T>
    void operator()(T& tcf) const
    {
	tcf(std::make_pair(sample.begin(), q_vector.begin()), std::make_pair(sample.end(), q_vector.end()), tcf.result[block].begin());
    }

private:
    unsigned int block;
    T1 const& sample;
    T2 const& q_vector;
};

template <typename T1, typename T2>
tcf_correlate_block<T1, T2> tcf_correlate_block_gen(unsigned int block, T1 const& sample, T2 const& q_vector)
{
    return tcf_correlate_block<T1, T2>(block, sample, q_vector);
}

/**
 * retrieve name of a correlation function
 */
class tcf_name : public boost::static_visitor<char const*>
{
public:
    template <typename T>
    char const* operator()(T const& tcf) const
    {
	return tcf.name();
    }
};

/**
 * create correlation function HDF5 dataset
 */
class tcf_create_dataset : public boost::static_visitor<>
{
public:
    tcf_create_dataset(H5::H5File& file, bool binary) : file(file), binary(binary) {}

    template <typename T>
    void operator()(T& tcf) const
    {
	H5::Group root(file.openGroup("/"));
	tcf.dataset = create_dataset(root, tcf.name(), tcf.result);
    }

    static H5::DataSet create_dataset(H5::Group const& node, char const* name, tcf_unary_result_type const& result)
    {
	// extensible dataspace for unary correlation function results
	hsize_t dim[3] = { 0, result.shape()[1], 5 };
	hsize_t max_dim[3] = { H5S_UNLIMITED, result.shape()[1], 5 };
	hsize_t chunk_dim[3] = { 1, result.shape()[1], 3 };
	H5::DataSpace ds(3, dim, max_dim);
	H5::DSetCreatPropList cparms;
	cparms.setChunk(3, chunk_dim);

	return node.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds, cparms);
    }

    static H5::DataSet create_dataset(H5::Group const& node, char const* name, tcf_binary_result_type const& result)
    {
	// extensible dataspace for binary correlation function results
	hsize_t dim[4] = { result.shape()[2], 0, result.shape()[1], 6 };
	hsize_t max_dim[4] = { result.shape()[2], H5S_UNLIMITED, result.shape()[1], 6 };
	hsize_t chunk_dim[4] = { result.shape()[2], 1, result.shape()[1], 4 };
	H5::DataSpace ds(4, dim, max_dim);
	H5::DSetCreatPropList cparms;
	cparms.setChunk(4, chunk_dim);

	return node.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds, cparms);
    }

private:
    H5::H5File& file;
    bool const binary;
};

/**
 * allocate correlation functions results
 */
class tcf_allocate_results : public boost::static_visitor<>
{
public:
    tcf_allocate_results(unsigned int block_count, unsigned int block_size, unsigned int q_values)
	: block_count(block_count), block_size(block_size), q_values(q_values) {}

    template <typename T>
    void operator()(T& tcf) const
    {
	resize(tcf.result);
    }

    void resize(tcf_unary_result_type& result) const
    {
	result.resize(boost::extents[block_count][block_size]);
    }

    void resize(tcf_binary_result_type& result) const
    {
	result.resize(boost::extents[block_count][block_size][q_values]);
    }

private:
    unsigned int block_count;
    unsigned int block_size;
    unsigned int q_values;
};

/**
 * write correlation function results to HDF5 file
 */
class tcf_write_results : public boost::static_visitor<>
{
public:
    typedef boost::multi_array<double, 2> block_time_type;
    typedef std::vector<double> q_value_vector;

public:
    tcf_write_results(block_time_type const& block_time, q_value_vector const& q_value, unsigned int max_blocks)
	: block_time(block_time), q_value(q_value), max_blocks(max_blocks) {}

    template <typename T>
    void operator()(T& tcf) const
    {
	write(tcf.dataset, tcf.result);
    }

    void write(H5::DataSet& dataset, tcf_unary_result_type const& result) const
    {
	// dataset dimensions
	boost::array<hsize_t, 3> dim = {{ max_blocks, result.shape()[1], 5 }};
	// memory buffer for results
	boost::multi_array<double, 3> data(dim);

	for (unsigned int j = 0; j < dim[0]; ++j) {
	    for (unsigned int k = 0; k < dim[1]; ++k) {
		// time interval
		data[j][k][0] = block_time[j][k];
		// mean average
		data[j][k][1] = result[j][k].mean();
		// standard error of mean
		data[j][k][2] = result[j][k].err();
		// variance
		data[j][k][3] = result[j][k].var();
		// count
		data[j][k][4] = result[j][k].count();
	    }
	}
	dataset.extend(dim.c_array());
	// write results to HDF5 file
	dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    }

    void write(H5::DataSet& dataset, tcf_binary_result_type const& result) const
    {
	// dataset dimensions
	boost::array<hsize_t, 4> dim = {{ result.shape()[2], max_blocks, result.shape()[1], 6 }};
	// memory buffer for results
	boost::multi_array<double, 4> data(dim);

	for (unsigned int j = 0; j < dim[0]; ++j) {
	    for (unsigned int k = 0; k < dim[1]; ++k) {
		for (unsigned int l = 0; l < dim[2]; ++l) {
		    // q-value
		    data[j][k][l][0] = q_value[j];
		    // time interval
		    data[j][k][l][1] = block_time[k][l];
		    // mean average
		    data[j][k][l][2] = result[k][l][j].mean();
		    // standard error of mean
		    data[j][k][l][3] = result[k][l][j].err();
		    // variance
		    data[j][k][l][4] = result[k][l][j].var();
		    // count
		    data[j][k][l][5] = result[k][l][j].count();
		}
	    }
	}
	dataset.extend(dim.c_array());
	// write results to HDF5 file
	dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    }

private:
    block_time_type const& block_time;
    q_value_vector const& q_value;
    unsigned int max_blocks;
};

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_TCF_VISITOR_HPP */
