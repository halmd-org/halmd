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

#ifndef MDSIM_TCF_HPP
#define MDSIM_TCF_HPP

#include <H5Cpp.h>
#include <boost/array.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/multi_array.hpp>
#include <boost/variant.hpp>
#include "accumulator.hpp"

namespace mdsim {

/** correlation function result types */
typedef boost::multi_array<accumulator<double>, 2> tcf_unary_result_type;
typedef boost::multi_array<accumulator<double>, 3> tcf_binary_result_type;

/**
 * mean-square displacement
 */
struct mean_square_displacement
{
    /** block sample results */
    tcf_unary_result_type result;
    /** HDF5 dataset */
    H5::DataSet dataset;

    char const* name() { return "MSD"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type::value_type::vector_type::const_iterator vector_const_iterator;
	typedef typename input_iterator::first_type::value_type::vector_type::value_type vector_type;
	typedef typename input_iterator::first_type::value_type::vector_type::value_type::value_type value_type;

	// iterate over phase space samples in block
	for (typename input_iterator::first_type it = first.first; it != last.first; ++it, ++result) {
	    // iterate over particle coordinates in current and first sample
	    for (vector_const_iterator r = it->r.begin(), r0 = first.first->r.begin(); r != it->r.end(); ++r, ++r0) {
		// displacement of particle
		vector_type dr = *r0 - *r;
		// accumulate square displacement
		*result += dr * dr;
	    }
	}
    }
};

/**
 * mean-quartic displacement
 */
struct mean_quartic_displacement
{
    /** block sample results */
    tcf_unary_result_type result;
    /** HDF5 dataset */
    H5::DataSet dataset;

    char const* name() { return "MQD"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type::value_type::vector_type::const_iterator vector_const_iterator;
	typedef typename input_iterator::first_type::value_type::vector_type::value_type vector_type;
	typedef typename input_iterator::first_type::value_type::vector_type::value_type::value_type value_type;

	// iterate over phase space samples in block
	for (typename input_iterator::first_type it = first.first; it != last.first; ++it, ++result) {
	    // iterate over particle coordinates in current and first sample
	    for (vector_const_iterator r = it->r.begin(), r0 = first.first->r.begin(); r != it->r.end(); ++r, ++r0) {
		// displacement of particle
		vector_type dr = *r0 - *r;
		// square displacement
		value_type rr = dr * dr;
		// accumulate quartic displacement
		*result += rr * rr;
	    }
	}
    }
};

/**
 * velocity autocorrelation
 */
struct velocity_autocorrelation
{
    /** block sample results */
    tcf_unary_result_type result;
    /** HDF5 dataset */
    H5::DataSet dataset;

    char const* name() { return "VAC"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type::value_type::vector_type::const_iterator vector_const_iterator;

	// iterate over phase space samples in block
	for (typename input_iterator::first_type it = first.first; it != last.first; ++it, ++result) {
	    // iterate over particle velocities in current and first sample
	    for (vector_const_iterator v = it->v.begin(), v0 = first.first->v.begin(); v != it->v.end(); ++v, ++v0) {
		// accumulate velocity autocorrelation
		*result += *v0 * *v;
	    }
	}
    }
};

/**
 * intermediate scattering function
 */
struct intermediate_scattering_function
{
    /** block sample results */
    tcf_binary_result_type result;
    /** HDF5 dataset */
    H5::DataSet dataset;

    char const* name() { return "ISF"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typename input_iterator::first_type::value_type::density_vector_type::const_iterator rho, rho0;
	typename output_iterator::value_type::iterator k;

	// iterate over phase space samples in block
	for (typename input_iterator::first_type it = first.first; it != last.first; ++it, ++result) {
	    // iterate over Fourier transformed densities in current and first sample
	    for (rho = it->rho.begin(), rho0 = first.first->rho.begin(), k = result->begin(); rho != it->rho.end(); ++rho, ++rho0, ++k) {
		// accumulate intermediate scattering function
		for (unsigned int d = 0; d < rho->first.size(); ++d) {
		    *k += rho->first[d] * rho0->first[d] + rho->second[d] * rho0->second[d];
		}
	    }
	}
    }
};

/**
 * self-intermediate scattering function
 */
struct self_intermediate_scattering_function
{
    /** block sample results */
    tcf_binary_result_type result;
    /** HDF5 dataset */
    H5::DataSet dataset;

    char const* name() { return "SISF"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typename input_iterator::first_type::value_type::vector_type::const_iterator r, r0;
	typename input_iterator::second_type q;
	typename input_iterator::first_type::value_type::vector_type::value_type f;
	unsigned int i;
	typename output_iterator::value_type::iterator k;

	// iterate over phase space samples in block
	for (typename input_iterator::first_type it = first.first; it != last.first; ++it, ++result) {
	    // iterate over q-values
	    for (q = first.second, k = result->begin(); q != last.second; ++q, ++k) {
		// iterate over particle positions in current and first sample
		for (f = 0, r = it->r.begin(), r0 = first.first->r.begin(), i = 0; r != it->r.end(); ++r, ++r0, ++i) {
		    f += (cos((*r - *r0) * *q) - f) / (i + 1);
		}
		// accumulate self-intermediate scattering function over mutually orthogonal q-vectors
		for (i = 0; i < f.size(); ++i) {
		    *k += f[i];
		}
	    }
	}
    }
};

/** correlation function types */
typedef boost::mpl::vector<mean_square_displacement, mean_quartic_displacement, velocity_autocorrelation, intermediate_scattering_function, self_intermediate_scattering_function> tcf_types;

/**
 * apply correlation function to block of phase space samples
 */
template <typename T1, typename T2>
class tcf_correlate_block : public boost::static_visitor<>
{
public:
    tcf_correlate_block(unsigned int const& block, T1 const& sample, T2 const& q_vector) : block(block), sample(sample), q_vector(q_vector) {}

    template <typename T>
    void operator()(T& tcf) const
    {
	tcf(std::make_pair(sample.begin(), q_vector.begin()), std::make_pair(sample.end(), q_vector.end()), tcf.result[block].begin());
    }

private:
    unsigned int const& block;
    T1 const& sample;
    T2 const& q_vector;
};

template <typename T1, typename T2>
tcf_correlate_block<T1, T2> tcf_correlate_block_gen(unsigned int const& block, T1 const& sample, T2 const& q_vector)
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
    char const* operator()(T& tcf) const
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
    tcf_create_dataset(H5::H5File& file) : file(file) {}

    template <typename T>
    void operator()(T& tcf) const
    {
	create_dataset(tcf.dataset, tcf.result, tcf.name());
    }

    void create_dataset(H5::DataSet& dataset, tcf_unary_result_type const& result, char const* name) const
    {
	// extensible dataspace for correlation function results
	hsize_t dim[3] = { 0, result.shape()[1], 5 };
	hsize_t max_dim[3] = { H5S_UNLIMITED, result.shape()[1], 5 };
	hsize_t chunk_dim[3] = { 1, result.shape()[1], 3 };
	H5::DataSpace ds(3, dim, max_dim);
	H5::DSetCreatPropList cparms;
	cparms.setChunk(3, chunk_dim);

	// create dataset
	dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds, cparms);
    }

    void create_dataset(H5::DataSet& dataset, tcf_binary_result_type const& result, char const* name) const
    {
	// extensible dataspace for binary correlation function results
	hsize_t dim[4] = { result.shape()[2], 0, result.shape()[1], 6 };
	hsize_t max_dim[4] = { result.shape()[2], H5S_UNLIMITED, result.shape()[1], 6 };
	hsize_t chunk_dim[4] = { result.shape()[2], 1, result.shape()[1], 4 };
	H5::DataSpace ds(4, dim, max_dim);
	H5::DSetCreatPropList cparms;
	cparms.setChunk(4, chunk_dim);

	// create dataset
	dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds, cparms);
    }
    
private:
    H5::H5File& file;
};

/**
 * allocate correlation functions results
 */
class tcf_allocate_results : public boost::static_visitor<>
{
public:
    tcf_allocate_results(unsigned int const& block_count, unsigned int const& block_size, unsigned int const& q_values) : block_count(block_count), block_size(block_size), q_values(q_values) {}

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
    unsigned int const& block_count;
    unsigned int const& block_size;
    unsigned int const& q_values;
};

/**
 * write correlation function results to HDF5 file
 */
class tcf_write_results : public boost::static_visitor<>
{
public:
    tcf_write_results(boost::multi_array<double, 2> const& block_time, std::vector<double> const& q_vector, unsigned int const& max_blocks) : block_time(block_time), q_vector(q_vector), max_blocks(max_blocks) {}

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
		    data[j][k][l][0] = q_vector[j];
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
    boost::multi_array<double, 2> const& block_time;
    std::vector<double> const& q_vector;
    unsigned int const& max_blocks;
};

} // namespace mdsim

#endif /* ! MDSIM_TCF_HPP */
