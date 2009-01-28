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

#ifndef LJGPU_SAMPLE_TCF_HOST_HPP
#define LJGPU_SAMPLE_TCF_HOST_HPP

#include <algorithm>
#include <boost/mpl/vector.hpp>
#include <boost/variant.hpp>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/mdsim/sample.hpp>
#include <ljgpu/sample/tcf_base.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <vector>

namespace ljgpu {

/**
 * Phase space sample for evaluating correlation functions
 */
template <int dimension>
struct tcf_host_sample : public tcf_sample<dimension>
{
    typedef tcf_sample<dimension> _Base;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::q_value_vector q_value_vector;
    typedef typename _Base::q_vector_vector q_vector_vector;
    typedef typename _Base::density_pair density_pair;
    typedef typename _Base::density_vector density_vector;
    typedef typename _Base::density_vector_vector density_vector_vector;
    typedef std::vector<vector_type> sample_vector;

    tcf_host_sample() {}
    tcf_host_sample(boost::shared_ptr<sample_vector> r, boost::shared_ptr<sample_vector> v): r(r), v(v) {}

    /**
     * initialise phase space sample
     */
    void operator()(q_vector_vector const& q)
    {
	typename q_vector_vector::const_iterator q0;
	typename q_vector_vector::value_type::const_iterator q1;
	typename density_vector_vector::iterator rho0;
	typename density_vector_vector::value_type::iterator rho1;
	typename sample_vector::const_iterator r0;

	// allocate memory for Fourier transformed densities
	rho = boost::shared_ptr<density_vector_vector>(new density_vector_vector(q.size()));
	for (q0 = q.begin(), rho0 = rho->begin(); q0 != q.end(); ++q0, ++rho0) {
	    rho0->assign(q0->size(), density_pair(0, 0));
	}
	// spatial Fourier transformation
	for (q0 = q.begin(), rho0 = rho->begin(); q0 != q.end(); ++q0, ++rho0) {
	    for (q1 = q0->begin(), rho1 = rho0->begin(); q1 != q0->end(); ++q1, ++rho1) {
		for (r0 = r->begin(); r0 != r->end(); ++r0) {
		    double const r_q = (*r0) * (*q1);
		    rho1->first += std::cos(r_q);
		    rho1->second += std::sin(r_q);
		}
	    }
	}
    }

    /** particle positions */
    boost::shared_ptr<sample_vector> r;
    /** particle velocities */
    boost::shared_ptr<sample_vector> v;
    /** Fourier transformed density for different |q| values and vectors */
    boost::shared_ptr<density_vector_vector> rho;
};

template <>
struct correlation_function<tcf_host_sample>
{
    /** HDF5 dataset */
    H5::DataSet dataset;
    /** particle type */
    size_t type;
};

/**
 * mean-square displacement
 */
template <>
struct mean_square_displacement<tcf_host_sample> : correlation_function<tcf_host_sample>
{
    /** block sample results */
    tcf_unary_result_type result;

    char const* name() const { return "MSD"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type::value_type sample_type;
	typedef typename sample_type::sample_vector::const_iterator vector_const_iterator;
	typedef typename sample_type::sample_vector::value_type vector_type;

	// iterate over phase space samples in block
	for (sample_iterator it = first.first; it != last.first; ++it, ++result) {
	    // iterate over particle coordinates in current and first sample
	    for (vector_const_iterator r = (*it)[type].r->begin(), r0 = (*first.first)[type].r->begin(); r != (*it)[type].r->end(); ++r, ++r0) {
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
template <>
struct mean_quartic_displacement<tcf_host_sample> : correlation_function<tcf_host_sample>
{
    /** block sample results */
    tcf_unary_result_type result;

    char const* name() const { return "MQD"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type::value_type sample_type;
	typedef typename sample_type::sample_vector::const_iterator vector_const_iterator;
	typedef typename sample_type::sample_vector::value_type vector_type;
	typedef typename sample_type::sample_vector::value_type::value_type value_type;

	// iterate over phase space samples in block
	for (sample_iterator it = first.first; it != last.first; ++it, ++result) {
	    // iterate over particle coordinates in current and first sample
	    for (vector_const_iterator r = (*it)[type].r->begin(), r0 = (*first.first)[type].r->begin(); r != (*it)[type].r->end(); ++r, ++r0) {
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
template <>
struct velocity_autocorrelation<tcf_host_sample> : correlation_function<tcf_host_sample>
{
    /** block sample results */
    tcf_unary_result_type result;

    char const* name() const { return "VAC"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type::value_type sample_type;
	typedef typename sample_type::sample_vector::const_iterator vector_const_iterator;

	// iterate over phase space samples in block
	for (sample_iterator it = first.first; it != last.first; ++it, ++result) {
	    // iterate over particle velocities in current and first sample
	    for (vector_const_iterator v = (*it)[type].v->begin(), v0 = (*first.first)[type].v->begin(); v != (*it)[type].v->end(); ++v, ++v0) {
		// accumulate velocity autocorrelation
		*result += *v0 * *v;
	    }
	}
    }
};

/**
 * intermediate scattering function
 */
template <>
struct intermediate_scattering_function<tcf_host_sample> : correlation_function<tcf_host_sample>
{
    /** block sample results */
    tcf_binary_result_type result;

    char const* name() const { return "ISF"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type::value_type sample_type;
	typedef typename sample_type::density_vector_vector::const_iterator density_vector_iterator;
	typedef typename sample_type::density_vector_vector::value_type::const_iterator density_iterator;
	typedef typename output_iterator::value_type::iterator q_value_result_iterator;

	for (sample_iterator it = first.first; it != last.first; ++it, ++result) {
	    q_value_result_iterator k = (*result).begin();
	    density_vector_iterator j0 = (*first.first)[type].rho->begin();
	    for (density_vector_iterator j = (*it)[type].rho->begin(); j != (*it)[type].rho->end(); ++j, ++j0, ++k) {
		density_iterator rho0 = (*j0).begin();
		for (density_iterator rho = (*j).begin(); rho != (*j).end(); ++rho, ++rho0) {
		    *k += (rho->first * rho0->first + rho->second * rho0->second) / (*it)[type].r->size();
		}
	    }
	}
    }
};

/**
 * self-intermediate scattering function
 */
template <>
struct self_intermediate_scattering_function<tcf_host_sample> : correlation_function<tcf_host_sample>
{
    /** block sample results */
    tcf_binary_result_type result;

    char const* name() const { return "SISF"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::first_type sample_iterator;
	typedef typename sample_iterator::value_type::value_type sample_type;
	typedef typename sample_type::sample_vector::const_iterator position_iterator;
	typedef typename sample_type::vector_type vector_type;
	typedef typename input_iterator::second_type q_value_iterator;
	typedef typename q_value_iterator::value_type::const_iterator q_vector_iterator;
	typedef typename output_iterator::value_type::iterator q_value_result_iterator;

	for (sample_iterator it = first.first; it != last.first; ++it, ++result) {
	    q_value_result_iterator k = (*result).begin();
	    for (q_value_iterator j = first.second; j != last.second; ++j, ++k) {
		for (q_vector_iterator q = (*j).begin(); q != (*j).end(); ++q) {
		    double value = 0;
		    size_t count = 0;
		    position_iterator r0 = (*first.first)[type].r->begin();
		    for (position_iterator r = (*it)[type].r->begin(); r != (*it)[type].r->end(); ++r, ++r0) {
			value += (std::cos((*r - *r0) * (*q)) - value) / (count + 1);
		    }
		    // result is normalised as we computed the average above
		    *k += value;
		}
	    }
	}
    }
};

/** correlation function types */
typedef boost::mpl::vector<mean_square_displacement<tcf_host_sample> > _tcf_host_types_0;
typedef boost::mpl::push_back<_tcf_host_types_0, mean_quartic_displacement<tcf_host_sample> >::type _tcf_host_types_1;
typedef boost::mpl::push_back<_tcf_host_types_1, velocity_autocorrelation<tcf_host_sample> >::type _tcf_host_types_2;
typedef boost::mpl::push_back<_tcf_host_types_2, intermediate_scattering_function<tcf_host_sample> >::type _tcf_host_types_3;
typedef boost::mpl::push_back<_tcf_host_types_3, self_intermediate_scattering_function<tcf_host_sample> >::type tcf_host_types;

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_TCF_HOST_HPP */
