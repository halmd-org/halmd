/*
 * Copyright © 2008-2014 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

// uncomment this line for a thread-divergent, slower implementation
// of the original thermostat introduced by H. C. Andersen (1978).
// #define USE_ORIGINAL_ANDERSEN_THERMOSTAT

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/integrators/verlet_nvt_andersen_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/random/gpu/normal_distribution.cuh>
#include <halmd/random/gpu/random_number_generator.cuh>
#include <halmd/utility/gpu/thread.cuh>

#if __CUDA_ARCH__ < 120
# define USE_ORIGINAL_ANDERSEN_THERMOSTAT
#endif

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {
namespace verlet_nvt_andersen_kernel {

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type, typename gpu_vector_type>
__global__ void integrate(
    float4* g_position
  , gpu_vector_type* g_image
  , float4* g_velocity
  , gpu_vector_type const* g_force
  , float timestep
  , fixed_vector<float, dimension> box_length
)
{
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef fixed_vector<float, dimension> float_vector_type;

    // kernel execution parameters
    unsigned int const thread = GTID;
    unsigned int const nthread = GTDIM;

    // read position, species, velocity, mass, image, force from global memory
    vector_type r, v;
    unsigned int species;
    float mass;
#ifdef USE_VERLET_DSFUN
    tie(r, species) <<= tie(g_position[thread], g_position[thread + nthread]);
    tie(v, mass) <<= tie(g_velocity[thread], g_velocity[thread + nthread]);
#else
    tie(r, species) <<= g_position[thread];
    tie(v, mass) <<= g_velocity[thread];
#endif
    float_vector_type f = g_force[thread];

    // advance position by full step, velocity by half step
    v += f * (timestep / 2) / mass;
    r += v * timestep;
    float_vector_type image = box_kernel::reduce_periodic(r, box_length);

    // store position, species, velocity, mass, image in global memory
#ifdef USE_VERLET_DSFUN
    tie(g_position[thread], g_position[thread + nthread]) <<= tie(r, species);
    tie(g_velocity[thread], g_velocity[thread + nthread]) <<= tie(v, mass);
#else
    g_position[thread] <<= tie(r, species);
    g_velocity[thread] <<= tie(v, mass);
#endif
    if (!(image == float_vector_type(0))) {
        g_image[thread] = image + static_cast<float_vector_type>(g_image[thread]);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 *
 * CUDA execution dimensions must agree with random number generator
 *
 * @param g_velocity particle velocities (array of size \code{} 2 * nplace \endcode for dsfloat arithmetic)
 * @param g_force particle forces (array of size \code{} nplace \endcode)
 * @param timestep integration time-step
 * @param sqrt_temperature square-root of heat bath temperature
 * @param coll_prob collision probability with heat bath
 * @param npart number of particles
 * @param nplace number of placeholder particles
 * @param rng random number generator
 */
template <int dimension, typename float_type, typename gpu_vector_type, typename rng_type>
__global__ void finalize(
    float4* g_velocity
  , gpu_vector_type const* g_force
  , float timestep
  , float sqrt_temperature
  , float coll_prob
  , unsigned int npart
  , unsigned int nplace
  , rng_type rng
)
{
    // read random number generator state from global device memory
    typename rng_type::state_type state = rng[GTID];

    // cache second normal variate for odd dimensions
    bool cached = false;
    float_type cache;

    // a heat bath collision is performed for all particles of a warp and
    // the collision probability thus requires a correction
    //
    // the probability to have no collision in _all_ threads of a warp
    // must equal the probability to have no collision in a single thread
    // pow(1 - q, warpSize) = 1 - coll_prob_
#ifndef USE_ORIGINAL_ANDERSEN_THERMOSTAT
    float q = 1 - __powf(1 - coll_prob, 1.f / warpSize);  //< powf() is significantly slower
#endif

    for (uint i = GTID; i < npart; i += GTDIM) {
        // read velocity, mass from global device memory
        fixed_vector<float_type, dimension> v;
        float mass;
#ifdef USE_VERLET_DSFUN
        tie(v, mass) <<= tie(g_velocity[i], g_velocity[i + nplace]);
#else
        tie(v, mass) <<= g_velocity[i];
#endif

        // is this a deterministic step?
        //
        // to avoid divergent threads within a warp, the stochastic coupling
        // is performed for all particles in the thread as soon as one thread requests it
#ifdef USE_ORIGINAL_ANDERSEN_THERMOSTAT
        if (uniform(rng, state) > coll_prob) {
#else
        if (__all(uniform(rng, state) > q)) {
#endif
            // read force from global device memory
            fixed_vector<float, dimension> f = g_force[i];
            // update velocity
            v += f * (timestep / 2) / mass;
        }
        // stochastic coupling with heat bath
        else {
            // parameters for normal distribution
            float const mean = 0;
            float const sigma = sqrt_temperature;

            // assign random velocity according to Maxwell-Boltzmann distribution
            for (uint j = 0; j < dimension - 1; j += 2) {
                tie(v[j], v[j + 1]) = normal(rng, state, mean, sigma);
            }
            if (dimension % 2) {
                if ((cached = !cached)) {
                    tie(v[dimension - 1], cache) = normal(rng, state, mean, sigma);
                }
                else {
                    v[dimension - 1] = cache;
                }
            }
        }

        // write velocity, mass to global device memory
#ifdef USE_VERLET_DSFUN
        tie(g_velocity[i], g_velocity[i + nplace]) <<= tie(v, mass);
#else
        g_velocity[i] <<= tie(v, mass);
#endif
    }

    // store random number generator state in global device memory
    rng[GTID] = state;
}

} // namespace verlet_nvt_andersen_kernel

template <int dimension, typename rng_type>
verlet_nvt_andersen_wrapper<dimension, rng_type> const
verlet_nvt_andersen_wrapper<dimension, rng_type>::kernel = {
#ifdef USE_VERLET_DSFUN
    verlet_nvt_andersen_kernel::integrate<dimension, dsfloat>
  , verlet_nvt_andersen_kernel::finalize<dimension, dsfloat>
#else
    verlet_nvt_andersen_kernel::integrate<dimension, float>
  , verlet_nvt_andersen_kernel::finalize<dimension, float>
#endif
};

template class verlet_nvt_andersen_wrapper<3, random::gpu::rand48_rng>;
template class verlet_nvt_andersen_wrapper<2, random::gpu::rand48_rng>;

} // namespace mdsim
} // namespace gpu
} // namespace integrators
} // namespace halmd
