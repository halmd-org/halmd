/* Trajectory correlation functions
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
 * generic correlation function type
 */
typedef boost::mpl::vector<mean_square_displacement> _tcf_types_0;
typedef boost::mpl::push_back<_tcf_types_0, mean_quartic_displacement>::type _tcf_types_1;
typedef boost::mpl::push_back<_tcf_types_1, velocity_autocorrelation>::type tcf_types;

/**
 * apply generic correlation function to block of phase space samples
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
