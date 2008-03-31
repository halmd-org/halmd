/* Common random number generator routines for CUDA device functions
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

#ifndef MDSIM_GPU_RNG_H
#define MDSIM_GPU_RNG_H

//
// Define the following macros before including this header:
//
//   __rng		: namespace of the generator [rand48]
//   __rng_state_t	: generator state type [ushort3]
//

namespace __rng
{

/**
 * returns uniform random number
 */
__inline__ __device__ float uniform(__rng_state_t& state);


/**
 * returns n-dimensional random unit vector
 */
template <typename T>
__inline__ __device__ T unit_vector(__rng_state_t& state);


/**
 * returns 2-dimensional random unit vector
 */
template <>
__inline__ __device__ float2 unit_vector(__rng_state_t& state)
{
    float r = 2. * M_PI * uniform(state);
    return make_float2(cosf(r), sinf(r));
}


/**
 * returns 3-dimensional random unit vector
 *
 * The following method requires an average of 8/Pi =~ 2.55
 * uniform random numbers. It is described in
 *
 * G. Marsaglia, Choosing a Point from the Surface of a Sphere,
 * The Annals of Mathematical Statistics, 1972, 43, 645-646
 *
 * http://projecteuclid.org/euclid.aoms/1177692644
 */
template <>
__inline__ __device__ float3 unit_vector(__rng_state_t& state)
{
    float3 v;
    float r;

    do {
	v.x = 2. * uniform(state) - 1.;
	v.y = 2. * uniform(state) - 1.;
	r = v.x * v.x + v.y * v.y;
    }
    while (r >= 1.);

    v.z = 1. - 2. * r;
    r = 2. * sqrtf(1. - r);
    v.x *= r;
    v.y *= r;
    return v;
}

} // namespace __rng

#endif /* ! MDSIM_GPU_RNG_H */
