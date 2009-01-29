/* Time correlation functions
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

#ifndef LJGPU_SAMPLE_TCF_HPP
#define LJGPU_SAMPLE_TCF_HPP

#include <H5Cpp.h>
#include <algorithm>
#include <boost/array.hpp>
#include <boost/assign.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/multi_array.hpp>
#include <boost/variant.hpp>
#include <ljgpu/math/accum.hpp>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/mdsim/sample.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <vector>

namespace ljgpu {

/**
 * Phase space sample for evaluating correlation functions
 */
template <int dimension>
struct uniform_tcf_sample
{
    typedef vector<double, dimension> vector_type;
    typedef std::vector<vector_type> sample_vector;
    typedef std::vector<double> q_value_vector;
    typedef std::vector<std::vector<vector_type> > q_vector_vector;

    /** real and imaginary components of Fourier transformed density rho(q) */
    typedef std::pair<double, double> density_vector_pair;
    /** vector of Fourier transformed densities for different q-values */
    typedef std::vector<std::vector<density_vector_pair> > density_vector_vector;

    /**
     * initialise phase space sample
     */
    template <typename position_sample_vector, typename velocity_sample_vector>
    void operator()(position_sample_vector const& _r, velocity_sample_vector const& _v, q_vector_vector const& q)
    {
	// copy particle positions
	r.assign(_r.begin(), _r.end());
	// copy particle velocities
	v.assign(_v.begin(), _v.end());

	rho.resize(q.size());
	for (size_t i = 0; i < rho.size(); ++i) {
	    rho[i].assign(q[i].size(), density_vector_pair(0, 0));
	}

	// spatial Fourier transformation
	double const norm = sqrt(r.size());
	for (size_t i = 0; i < q.size(); ++i) {
	    for (size_t j = 0; j < q[i].size(); ++j) {
		for (size_t k = 0; k < r.size(); ++k) {
		    double const value = r[k] * q[i][j];
		    // sum over differences to maintain accuracy with small and large values
		    rho[i][j].first += (std::cos(value) - rho[i][j].first) / (k + 1);
		    rho[i][j].second += (std::sin(value) - rho[i][j].second) / (k + 1);
		}
		// multiply with norm as we compute average value above
		rho[i][j].first *= norm;
		rho[i][j].second *= norm;
	    }
	}
    }

    /** particle positions */
    sample_vector r;
    /** particle velocities */
    sample_vector v;
    /** spatially Fourier transformed density for given q-values */
    density_vector_vector rho;
};

template <int dimension>
struct tcf_sample : boost::array<uniform_tcf_sample<dimension>, 2>
{
    typedef boost::array<uniform_tcf_sample<dimension>, 2> _Base;
    typedef typename _Base::value_type::vector_type vector_type;
    typedef typename _Base::value_type::sample_vector sample_vector;
    typedef typename _Base::value_type::q_value_vector q_value_vector;
    typedef typename _Base::value_type::q_vector_vector q_vector_vector;
    typedef typename _Base::value_type::density_vector_pair density_vector_pair;
    typedef typename _Base::value_type::density_vector_vector density_vector_vector;

    template <typename sample_type>
    void operator()(sample_type const& sample, q_vector_vector const& q)
    {
	for (size_t i = 0; i < sample.size(); ++i) {
	    (*this)[i](sample[i].r, sample[i].v, q);
	}
    }
};

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
    /** particle types */
    particle_type a, b;

    mean_square_displacement() : a(PART_A), b(PART_A) {}
    mean_square_displacement(particle_type a, particle_type b) : a(a), b(b) {}

    char const* name() const { return "MSD"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type::value_type::sample_vector::const_iterator vector_const_iterator;
	typedef typename input_iterator::first_type::value_type::sample_vector::value_type vector_type;

	// iterate over phase space samples in block
	for (typename input_iterator::first_type it = first.first; it != last.first; ++it, ++result) {
	    // iterate over particle coordinates in current and first sample
	    for (vector_const_iterator r = (*it)[b].r.begin(), r0 = (*first.first)[a].r.begin(); r != (*it)[b].r.end(); ++r, ++r0) {
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
    /** particle types */
    particle_type a, b;

    mean_quartic_displacement() : a(PART_A), b(PART_A) {}
    mean_quartic_displacement(particle_type a, particle_type b) : a(a), b(b) {}

    char const* name() const { return "MQD"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type::value_type::sample_vector::const_iterator vector_const_iterator;
	typedef typename input_iterator::first_type::value_type::sample_vector::value_type vector_type;
	typedef typename input_iterator::first_type::value_type::sample_vector::value_type::value_type value_type;

	// iterate over phase space samples in block
	for (typename input_iterator::first_type it = first.first; it != last.first; ++it, ++result) {
	    // iterate over particle coordinates in current and first sample
	    for (vector_const_iterator r = (*it)[b].r.begin(), r0 = (*first.first)[a].r.begin(); r != (*it)[b].r.end(); ++r, ++r0) {
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
    /** particle types */
    particle_type a, b;

    velocity_autocorrelation() : a(PART_A), b(PART_A) {}
    velocity_autocorrelation(particle_type a, particle_type b) : a(a), b(b) {}

    char const* name() const { return "VAC"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type::value_type::sample_vector::const_iterator vector_const_iterator;

	// iterate over phase space samples in block
	for (typename input_iterator::first_type it = first.first; it != last.first; ++it, ++result) {
	    // iterate over particle velocities in current and first sample
	    for (vector_const_iterator v = (*it)[b].v.begin(), v0 = (*first.first)[a].v.begin(); v != (*it)[b].v.end(); ++v, ++v0) {
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
    /** particle types */
    particle_type a, b;

    intermediate_scattering_function() : a(PART_A), b(PART_A) {}
    intermediate_scattering_function(particle_type a, particle_type b) : a(a), b(b) {}

    char const* name() const { return "ISF"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type sample_type;
	typedef typename sample_type::density_vector_vector::const_iterator density_vector_iterator;
	typedef typename sample_type::density_vector_vector::value_type::const_iterator density_iterator;
	typedef typename output_iterator::value_type::iterator q_value_result_iterator;

	for (sample_iterator i = first.first; i != last.first; ++i, ++result) {
	    q_value_result_iterator k = (*result).begin();
	    density_vector_iterator j0 = (*first.first)[a].rho.begin();
	    for (density_vector_iterator j = (*i)[b].rho.begin(); j != (*i)[b].rho.end(); ++j, ++j0, ++k) {
		density_iterator rho0 = (*j0).begin();
		for (density_iterator rho = (*j).begin(); rho != (*j).end(); ++rho, ++rho0, ++k) {
		    *k += rho->first * rho0->first + rho->second * rho0->second;
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
    /** particle types */
    particle_type a, b;

    self_intermediate_scattering_function() : a(PART_A), b(PART_A) {}
    self_intermediate_scattering_function(particle_type a, particle_type b) : a(a), b(b) {}

    char const* name() const { return "SISF"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type sample_type;
	typedef typename sample_type::sample_vector::const_iterator position_iterator;
	typedef typename sample_type::vector_type vector_type;
	typedef typename input_iterator::second_type q_value_iterator;
	typedef typename q_value_iterator::value_type::const_iterator q_vector_iterator;
	typedef typename output_iterator::value_type::iterator q_value_result_iterator;
	typedef typename output_iterator::value_type::value_type::value_type q_value_result_value;

	for (sample_iterator i = first.first; i != last.first; ++i, ++result) {
	    q_value_result_iterator k = (*result).begin();
	    for (q_value_iterator j = first.second; j != last.second; ++j, ++k) {
		for (q_vector_iterator q = (*j).begin(); q != (*j).end(); ++q) {
		    q_value_result_value value = 0;
		    position_iterator r0 = (*first.first)[a].r.begin();
		    for (position_iterator r = (*i)[b].r.begin(); r != (*i)[b].r.end(); ++r, ++r0) {
			value += std::cos((*r - *r0) * (*q));
		    }
		    *k += value / (*i)[b].r.size();
		}
	    }
	}
    }
};

/** correlation function types */
typedef boost::mpl::vector<mean_square_displacement> _tcf_types_0;
typedef boost::mpl::push_back<_tcf_types_0, mean_quartic_displacement>::type _tcf_types_1;
typedef boost::mpl::push_back<_tcf_types_1, velocity_autocorrelation>::type _tcf_types_2;
typedef boost::mpl::push_back<_tcf_types_2, intermediate_scattering_function>::type _tcf_types_3;
typedef boost::mpl::push_back<_tcf_types_3, self_intermediate_scattering_function>::type tcf_types;

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
	using namespace boost::assign;
	boost::array<char const*, 3> const str = list_of("AA")("AB")("BB");

	H5::Group root(file.openGroup("/"));
	if (binary) {
	    try {
		H5XX_NO_AUTO_PRINT(H5::GroupIException);
		root = root.openGroup(str[tcf.a + tcf.b]);
	    }
	    catch (H5::GroupIException const&) {
		root = root.createGroup(str[tcf.a + tcf.b]);
	    }
	}
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

#endif /* ! LJGPU_SAMPLE_TCF_HPP */
