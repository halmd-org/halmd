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
// requires boost 1.37.0 or patch from http://svn.boost.org/trac/boost/ticket/1852
#include <boost/circular_buffer.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/multi_array.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/variant.hpp>
#ifdef WITH_CUDA
# include <ljgpu/sample/tcf_gpu.hpp>
#endif
#include <ljgpu/sample/tcf_host.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <vector>

namespace ljgpu {

template <typename mdsim_backend>
class _tcf_sample_phase_space : public boost::static_visitor<>
{
public:
    _tcf_sample_phase_space(mdsim_backend const& fluid) : fluid(fluid) {}

    template <int dimension>
    void operator() (tcf_host_sample<dimension>& sample) const
    {
	typedef tcf_host_sample<dimension> sample_type;
	typedef typename sample_type::sample_vector sample_vector;
	typedef boost::shared_ptr<sample_vector> sample_ptr;
	sample.r = sample_ptr(new sample_vector(fluid.sample()[0].r.begin(), fluid.sample()[0].r.end()));
	sample.v = sample_ptr(new sample_vector(fluid.sample()[0].v.begin(), fluid.sample()[0].v.end()));
    }

#ifdef WITH_CUDA

    /**
     * sample from global device memory
     */
    template <int dimension>
    typename boost::enable_if<boost::is_base_of<ljfluid_impl_gpu_neighbour<dimension>, typename mdsim_backend::impl_type>, void>::type
    operator() (tcf_gpu_sample<dimension>& sample) const
    {
	trajectory_gpu_sample<dimension> sample_;
	fluid.sample(sample_);
	sample.r = sample_.r[0];
	sample.v = sample_.v[0];
    }

    /**
     * sample from host memory
     */
    template <int dimension>
    typename boost::enable_if<boost::mpl::not_<boost::is_base_of<ljfluid_impl_gpu_neighbour<dimension>, typename mdsim_backend::impl_type> >, void>::type
    operator() (tcf_gpu_sample<dimension>& sample) const
    {
	typedef tcf_gpu_sample<dimension> sample_type;
	typedef typename sample_type::gpu_sample_vector gpu_sample_vector;
	typedef boost::shared_ptr<gpu_sample_vector> gpu_sample_ptr;
	typedef typename sample_type::gpu_vector_type gpu_vector_type;

	cuda::host::vector<gpu_vector_type> r(fluid.sample()[0].r.size());
	cuda::host::vector<gpu_vector_type> v(fluid.sample()[0].v.size());
	std::copy(fluid.sample()[0].r.begin(), fluid.sample()[0].r.end(), r.begin());
	std::copy(fluid.sample()[0].v.begin(), fluid.sample()[0].v.end(), v.begin());
	sample.r = gpu_sample_ptr(new gpu_sample_vector(r.size()));
	sample.v = gpu_sample_ptr(new gpu_sample_vector(v.size()));
	cuda::copy(r, *sample.r);
	cuda::copy(v, *sample.v);
    }

#endif /* WITH_CUDA */

private:
    mdsim_backend const& fluid;
};

template <typename mdsim_backend>
_tcf_sample_phase_space<mdsim_backend> tcf_sample_phase_space(mdsim_backend const& fluid)
{
    return _tcf_sample_phase_space<mdsim_backend>(fluid);
}

template <typename U>
class _tcf_fourier_transform_sample : public boost::static_visitor<>
{
public:
    _tcf_fourier_transform_sample(U const& q_vector) : q_vector(q_vector) {}

    template <typename T>
    void operator() (T& sample) const
    {
	sample(q_vector);
    }

private:
    U const& q_vector;
};

template <typename U>
_tcf_fourier_transform_sample<U> tcf_fourier_transform_sample(U const& q_vector)
{
   return _tcf_fourier_transform_sample<U>(q_vector);
}

class tcf_block_add_sample : public boost::static_visitor<>
{
public:
    template <typename T, typename U>
    void operator() (T const&, U const&) const
    {
	throw std::runtime_error("block sample mismatch");
    }

    template <typename T>
    void operator() (boost::circular_buffer<T>& block, T const& sample) const
    {
	block.push_back(sample);
    }
};

class tcf_block_is_full : public boost::static_visitor<bool>
{
public:
    template <typename T>
    bool operator() (boost::circular_buffer<T> const& block) const
    {
	return block.full();
    }
};

class tcf_block_clear : public boost::static_visitor<>
{
public:
    template <typename T>
    void operator() (boost::circular_buffer<T>& block) const
    {
	block.clear();
    }
};

class tcf_block_pop_front : public boost::static_visitor<>
{
public:
    template <typename T>
    void operator() (boost::circular_buffer<T>& block) const
    {
	block.pop_front();
    }
};

class tcf_block_size : public boost::static_visitor<size_t>
{
public:
    template <typename T>
    size_t operator() (boost::circular_buffer<T> const& block) const
    {
	return block.size();
    }
};

/**
 * apply correlation function to block of phase space samples
 */
template <typename V>
class _tcf_correlate_block : public boost::static_visitor<>
{
public:
    _tcf_correlate_block(unsigned int block, V const& q_vector) : block(block), q_vector(q_vector) {}

    template <typename T, typename U>
    void operator()(T&, U&) const
    {
	throw std::runtime_error("correlation function mismatch");
    }

    template <typename T, template <int> class sample_type, int dimension>
    typename boost::enable_if<boost::is_base_of<correlation_function<sample_type>, T>, void>::type
    operator()(T& tcf, boost::circular_buffer<sample_type<dimension> >& sample) const
    {
	tcf(std::make_pair(sample.begin(), q_vector.begin()), std::make_pair(sample.end(), q_vector.end()), tcf.result[block].begin());
    }

private:
    unsigned int block;
    V const& q_vector;
};

template <typename T>
_tcf_correlate_block<T> tcf_correlate_block(unsigned int block, T const& q_vector)
{
    return _tcf_correlate_block<T>(block, q_vector);
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
