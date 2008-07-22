/* Phase space sample
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

#ifndef MDSIM_SAMPLE_HPP
#define MDSIM_SAMPLE_HPP

#include <algorithm>
#include <cmath>
#include <cuda_wrapper.hpp>
#include <vector>

namespace mdsim {

/**
 * MD simulation sample
 */
template <typename T>
struct trajectory_sample
{
    /** periodically reduced particle positions */
    std::vector<T> r;
    /** periodically extended particle positions */
    std::vector<T> R;
    /** particle velocities */
    std::vector<T> v;
    /** potential energy per particle */
    float en_pot;
    /** virial equation sum per particle */
    float virial;
};

/**
 * Phase space sample
 */
template <typename T>
struct phase_space_point
{
    // swappable host memory vector type
    typedef std::vector<T> vector_type;
    typedef typename T::value_type value_type;
    typedef typename std::vector<std::pair<T, T> > density_vector_type;

    phase_space_point(vector_type const& r, vector_type const& v, std::vector<value_type> const& q) : r(r), v(v), rho(q.size(), std::pair<T, T>(0, 0))
    {
	// spatial Fourier transformation
	for (size_t i = 0; i < r.size(); ++i) {
	    for (unsigned int j = 0; j < q.size(); ++j) {
		// compute averages to maintain accuracy with single precision floating-point
		rho[j].first += (cos(r[i] * q[j]) - rho[j].first) / (i + 1);
		rho[j].second += (sin(r[i] * q[j]) - rho[j].second) / (i + 1);
	    }
	}
	// normalize Fourier transformed density with N^(-1/2)
	const value_type n = std::sqrt(r.size());
	for (unsigned int j = 0; j < q.size(); ++j) {
	    // multiply averages with N^(+1/2)
	    rho[j].first *= n;
	    rho[j].second *= n;
	}
    }

    /** particle positions */
    vector_type r;
    /** particle velocities */
    vector_type v;
    /** spatially Fourier transformed density for given q-values */
    density_vector_type rho;
};

} // namespace mdsim

#endif /* ! MDSIM_SAMPLE_HPP */
