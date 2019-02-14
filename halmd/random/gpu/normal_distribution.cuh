/*
 * Copyright © 2007-2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_RANDOM_GPU_NORMAL_DISTRIBUTION_CUH
#define HALMD_RANDOM_GPU_NORMAL_DISTRIBUTION_CUH

#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace random {
namespace gpu {

//
// The Box-Muller transformation for generating random numbers
// in the normal distribution was originally described in
//
// G.E.P. Box and M.E. Muller, A Note on the Generation of
// Random Normal Deviates, The Annals of Mathematical Statistics,
// 1958, 29, p. 610-611
//
// Here, we use instead the faster polar method of the Box-Muller
// transformation, see
//
// D.E. Knuth, Art of Computer Programming, Volume 2: Seminumerical
// Algorithms, 3rd Edition, 1997, Addison-Wesley, p. 122
//

/**
 * generate pair of random numbers from Gaussian distribution
 */
template <typename RandomNumberGenerator>
inline __device__ tuple<float, float> normal(
    RandomNumberGenerator const& rng
  , typename RandomNumberGenerator::state_type& state
  , float mean = 0.0    // mean of distribution
  , float sigma = 1.0   // standard deviation
)
{
    float variate1, variate2, s;
    do {
        variate1 = 2 * uniform(rng, state) - 1;
        variate2 = 2 * uniform(rng, state) - 1;
        s = variate1 * variate1 + variate2 * variate2;
    } while (s >= 1);

    s = sigma * sqrtf(-2 * logf(s) / s);
    variate1 = s * variate1 + mean;
    variate2 = s * variate2 + mean;
    return make_tuple(variate1, variate2);
}

} // namespace random
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_NORMAL_DISTRIBUTION_CUH */
