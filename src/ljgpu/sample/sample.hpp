/* Phase space sample
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

#ifndef LJGPU_SAMPLE_SAMPLE_HPP
#define LJGPU_SAMPLE_SAMPLE_HPP

#include <algorithm>
#include <boost/function.hpp>
#include <vector>

namespace ljgpu {

/**
 * MD simulation sample
 */
template <typename vector_type>
struct trajectory_gpu_sample
{
    typedef std::vector<vector_type> sample_vector;
    typedef boost::function<void (sample_vector&, sample_vector&)> sample_visitor;

    /** periodically extended and reduced particle positions */
    // FIXME std::vector<std::pair<vector_type, vector_type> > r;
    sample_vector R;
    sample_vector r;
    /** particle velocities */
    sample_vector v;
    /** potential energy per particle */
    double en_pot;
    /** virial equation sum per particle */
    double virial;
};

template <typename vector_type>
struct trajectory_host_sample
{
    typedef std::vector<vector_type> sample_vector;
    typedef boost::function<void (sample_vector&, sample_vector&)> sample_visitor;

    /** periodically extended particle positions */
    sample_vector R;
    /** particle velocities */
    sample_vector v;
    /** potential energy per particle */
    double en_pot;
    /** virial equation sum per particle */
    double virial;
};

/**
 * Phase space sample for evaluating correlation functions
 */
template <typename vector_type>
struct correlation_sample
{
    // value types of these two vectors must match
    typedef std::vector<vector_type> sample_vector;
    typedef std::vector<typename vector_type::value_type> q_value_vector;
    // summing over large particle numbers requires double precision!
    typedef double density_value;
    /** real or imaginary component vector */
    typedef vector<density_value, vector_type::static_size> density_vector;
    /** real and imaginary components of Fourier transformed density rho(q) */
    typedef std::pair<density_vector, density_vector> density_vector_pair;
    /** vector of Fourier transformed densities for different q-values */
    typedef std::vector<density_vector_pair> density_vector_vector;

    /** particle positions */
    sample_vector r;
    /** particle velocities */
    sample_vector v;
    /** spatially Fourier transformed density for given q-values */
    density_vector_vector rho;

    /**
     * initialise phase space sample
     */
    correlation_sample(sample_vector const& r, sample_vector const& v, q_value_vector const& q)
    : r(r), v(v), rho(q.size(), density_vector_pair(0, 0))
    {
	// spatial Fourier transformation
	for (size_t i = 0; i < r.size(); ++i) {
	    for (size_t j = 0; j < q.size(); ++j) {
		// compute averages to maintain accuracy summing over small and large values
		density_vector r_q(r[i] * q[j]);
		rho[j].first += (cos(r_q) - rho[j].first) / (i + 1);
		rho[j].second += (sin(r_q) - rho[j].second) / (i + 1);
	    }
	}
	// normalize Fourier transformed density with N^(-1/2)
	density_value n = sqrt(r.size());
	for (size_t j = 0; j < q.size(); ++j) {
	    // therefore multiply averages with N^(+1/2)
	    rho[j].first *= n;
	    rho[j].second *= n;
	}
    }
};

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_SAMPLE_HPP */
