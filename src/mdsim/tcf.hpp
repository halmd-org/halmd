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

#include <boost/variant.hpp>
#include <boost/mpl/vector.hpp>


namespace mdsim {

/**
 * mean-square displacement
 */
struct mean_square_displacement
{
    char const* name() { return "MSD"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::value_type::vector_type::const_iterator vector_const_iterator;
	typedef typename input_iterator::value_type::vector_type::value_type vector_type;
	typedef typename input_iterator::value_type::vector_type::value_type::value_type value_type;

	// iterate over phase space samples in block, skipping first sample
	for (input_iterator it = first; ++it != last; ++result) {
	    // iterate over particle coordinates in current and first sample
	    for (vector_const_iterator r = it->r.begin(), r0 = first->r.begin(); r != it->r.end(); ++r, ++r0) {
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
    char const* name() { return "MQD"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::value_type::vector_type::const_iterator vector_const_iterator;
	typedef typename input_iterator::value_type::vector_type::value_type vector_type;
	typedef typename input_iterator::value_type::vector_type::value_type::value_type value_type;

	// iterate over phase space samples in block, skipping first sample
	for (input_iterator it = first; ++it != last; ++result) {
	    // iterate over particle coordinates in current and first sample
	    for (vector_const_iterator r = it->r.begin(), r0 = first->r.begin(); r != it->r.end(); ++r, ++r0) {
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
    char const* name() { return "VAC"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typedef typename input_iterator::value_type::vector_type::const_iterator vector_const_iterator;

	// iterate over phase space samples in block, skipping first sample
	for (input_iterator it = first; ++it != last; ++result) {
	    // iterate over particle velocities in current and first sample
	    for (vector_const_iterator v = it->v.begin(), v0 = first->v.begin(); v != it->v.end(); ++v, ++v0) {
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
    char const* name() { return "ISF"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typename input_iterator::value_type::density_vector_type::const_iterator rho, rho0;
	typename output_iterator::value_type::iterator k;

	// iterate over phase space samples in block
	for (input_iterator it = first; it != last; ++it, ++result) {
	    // iterate over Fourier transformed densities in current and first sample
	    for (rho = it->rho.begin(), rho0 = first->rho.begin(), k = result->begin(); rho != it->rho.end(); ++rho, ++rho0, ++k) {
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
    char const* name() { return "SISF"; }

    template <typename input_iterator, typename output_iterator>
    void operator()(input_iterator const& first, input_iterator const& last, output_iterator result)
    {
	typename input_iterator::value_type::vector_type::const_iterator r, r0;
	typename input_iterator::value_type::q_values_type::const_iterator q;
	typename output_iterator::value_type::iterator k;

	// iterate over phase space samples in block
	for (input_iterator it = first; it != last; ++it, ++result) {
	    // iterate over q-values
	    for (q = it->q.first, k = result->begin(); q != it->q.second; ++q, ++k) {
		typename input_iterator::value_type::vector_type::value_type sum(0);
		unsigned int i = 0;
		// iterate over particle positions in current and first sample
		for (r = it->r.begin(), r0 = first->r.begin(); r != it->r.end(); ++r, ++r0, ++i) {
		    sum += (cos((*r - *r0) * *q) - sum) / (i + 1);
		}
		// accumulate self-intermediate scattering function
		for (unsigned int d = 0; d < it->r.begin()->size(); ++d) {
		    *k += sum[d];
		}
	    }
	}
    }
};

/** correlation function type */
typedef boost::mpl::vector<mean_square_displacement> _tcf_types_0;
typedef boost::mpl::push_back<_tcf_types_0, mean_quartic_displacement>::type _tcf_types_1;
typedef boost::mpl::push_back<_tcf_types_1, velocity_autocorrelation>::type tcf_types;
/** binary correlation function type */
typedef boost::mpl::vector<intermediate_scattering_function> _qtcf_types_0;
typedef boost::mpl::push_back<_qtcf_types_0, self_intermediate_scattering_function>::type qtcf_types;

/**
 * apply correlation function to block of phase space samples
 */
template <typename T1, typename T2>
class tcf_apply_visitor : public boost::static_visitor<>
{
public:
    tcf_apply_visitor(T1 const& first, T1 const& last, T2 const& result) : first(first), last(last), result(result) { }

    template <typename T>
    void operator()(T& tcf) const
    {
	tcf(first, last, result);
    }

private:
    T1 const& first, last;
    T2 const& result;
};

template <typename T1, typename T2>
tcf_apply_visitor<T1, T2> tcf_apply_visitor_gen(T1 const& first, T1 const& last, T2 const& result)
{
    return tcf_apply_visitor<T1, T2>(first, last, result);
}

/**
 * retrieve name of a generic correlation function
 */
struct tcf_name_visitor : public boost::static_visitor<char const*>
{
    template <typename T>
    char const* operator()(T& tcf) const
    {
	return tcf.name();
    }
};

} // namespace mdsim

#endif /* ! MDSIM_TCF_HPP */
