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

#ifndef LJGPU_SAMPLE_TCF_BASE_HPP
#define LJGPU_SAMPLE_TCF_BASE_HPP

#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <ljgpu/math/accum.hpp>
#include <vector>

namespace ljgpu {

template <int dimension>
struct tcf_sample
{
    /** q vector in momentum space */
    typedef vector<double, dimension> vector_type;
    /** |q| value vector */
    typedef std::vector<double> q_value_vector;
    /** vector of q vectors for different |q| values */
    typedef std::vector<std::vector<vector_type> > q_vector_vector;
    /** real and imaginary components of Fourier transformed density rho(q) */
    typedef std::pair<double, double> density_pair;
    /** vector of Fourier transformed densities for |q| value */
    typedef std::vector<density_pair> density_vector;
    /** vector of density vectors */
    typedef std::vector<density_vector> density_vector_vector;
};

/** correlation function result types */
typedef boost::multi_array<accumulator<double>, 2> tcf_unary_result_type;
typedef boost::multi_array<accumulator<double>, 3> tcf_binary_result_type;

template <template <int> class sample_type>
struct correlation_function;

template <template <int> class sample_type>
struct mean_square_displacement;

template <template <int> class sample_type>
struct mean_quartic_displacement;

template <template <int> class sample_type>
struct velocity_autocorrelation;

template <template <int> class sample_type>
struct intermediate_scattering_function;

template <template <int> class sample_type>
struct self_intermediate_scattering_function;

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_TCF_BASE_HPP */
