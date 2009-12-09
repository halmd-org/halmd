/* Time correlation functions
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#ifndef HALMD_SAMPLE_TCF_HOST_HPP
#define HALMD_SAMPLE_TCF_HOST_HPP

#include <algorithm>
#include <boost/mpl/vector.hpp>
#include <boost/variant.hpp>
#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>
#include <halmd/mdsim/sample.hpp>
#include <halmd/sample/tcf_base.hpp>
#include <halmd/util/H5xx.hpp>
#include <vector>

namespace halmd {

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
    typedef typename _Base::isf_vector_vector isf_vector_vector;
    typedef typename _Base::density_pair density_pair;
    typedef typename _Base::density_vector density_vector;
    typedef typename _Base::density_vector_vector density_vector_vector;
    typedef typename _Base::virial_tensor virial_tensor;
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
        typename isf_vector_vector::iterator isf0;
        typename sample_vector::const_iterator r0;

        // allocate memory for Fourier transformed densities
        rho = boost::shared_ptr<density_vector_vector>(new density_vector_vector(q.size()));
        isf = boost::shared_ptr<isf_vector_vector>(new isf_vector_vector(q.size()));
        for (q0 = q.begin(), rho0 = rho->begin(), isf0 = isf->begin(); q0 != q.end(); ++q0, ++rho0, ++isf0) {
            rho0->assign(q0->size(), density_pair(0, 0));
            isf0->resize(q0->size());
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
    /** self-intermediate scattering function for different |q| values and vectors */
    boost::shared_ptr<isf_vector_vector> isf;
    /** off-diagonal elements of virial stress tensor */
    boost::shared_ptr<virial_tensor> virial;
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
        typedef typename sample_type::isf_vector_vector isf_vector_vector;
        typedef typename isf_vector_vector::iterator isf_vector_iterator;
        typedef typename isf_vector_vector::value_type::iterator isf_value_iterator;
        typedef typename input_iterator::second_type q_value_iterator;
        typedef typename q_value_iterator::value_type::const_iterator q_vector_iterator;
        typedef typename output_iterator::value_type::iterator q_value_result_iterator;
        typedef typename output_iterator::value_type::value_type::value_type q_value_result_value;

        for (sample_iterator it = first.first; it != last.first; ++it, ++result) {
            q_value_result_iterator k = (*result).begin();
            isf_vector_iterator isf0 = (*it)[type].isf->begin();
            for (q_value_iterator j = first.second; j != last.second; ++j, ++k, ++isf0) {
                isf_value_iterator isf = isf0->begin();
                for (q_vector_iterator q = (*j).begin(); q != (*j).end(); ++q, ++isf) {
                    q_value_result_value value = 0;
                    position_iterator r0 = (*first.first)[type].r->begin();
                    for (position_iterator r = (*it)[type].r->begin(); r != (*it)[type].r->end(); ++r, ++r0) {
                        value += std::cos((*r - *r0) * (*q));
                    }
                    *isf = value / (*it)[type].r->size();
                    *k += *isf;
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
typedef boost::mpl::push_back<_tcf_host_types_3, self_intermediate_scattering_function<tcf_host_sample> >::type _tcf_host_types_4;
typedef boost::mpl::push_back<_tcf_host_types_4, squared_self_intermediate_scattering_function<tcf_host_sample> >::type _tcf_host_types_5;
typedef boost::mpl::push_back<_tcf_host_types_5, virial_stress<tcf_host_sample> >::type tcf_host_types;

} // namespace halmd

#endif /* ! HALMD_SAMPLE_TCF_HOST_HPP */
