/*
 * Copyright © 2008-2010  Peter Colberg
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

// uncomment this line for a thread-divergent, slower implementation
// of the original thermostat introduced by H. C. Andersen (1978).
// #define USE_ORIGINAL_ANDERSEN_THERMOSTAT

#include <boost/mpl/if.hpp>

#include <halmd/mdsim/gpu/integrators/verlet_kernel.cuh>
#include <halmd/mdsim/gpu/integrators/verlet_nvt_andersen_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/random/gpu/normal_distribution.cuh>
#include <halmd/random/gpu/random_number_generator.cuh>
#include <halmd/utility/gpu/thread.cuh>

#if __CUDA_ARCH__ < 120
# define USE_ORIGINAL_ANDERSEN_THERMOSTAT
#endif

using namespace boost::mpl;
using namespace halmd::mdsim::gpu::particle_kernel;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {
namespace verlet_nvt_andersen_kernel {

/** integration time-step */
static __constant__ float timestep_;
/** square-root of heat bath temperature */
static __constant__ float sqrt_temperature_;
/** collision probability with heat bath */
static __constant__ float coll_prob_;

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <
    typename vector_type
  , typename vector_type_
  , typename gpu_vector_type
>
__global__ void _integrate(
    float4* g_r
  , gpu_vector_type* g_image
  , float4* g_v
  , gpu_vector_type const* g_f
  , vector_type_ box_length
)
{
    unsigned int const i = GTID;
    unsigned int const threads = GTDIM;
    unsigned int type, tag;
    vector_type r, v;
#ifdef USE_VERLET_DSFUN
    tie(r, type) = untagged<vector_type>(g_r[i], g_r[i + threads]);
    tie(v, tag) = untagged<vector_type>(g_v[i], g_v[i + threads]);
#else
    tie(r, type) = untagged<vector_type>(g_r[i]);
    tie(v, tag) = untagged<vector_type>(g_v[i]);
#endif
    vector_type_ image = g_image[i];
    vector_type_ f = g_f[i];

    verlet_kernel::integrate(r, image, v, f, timestep_, box_length);

#ifdef USE_VERLET_DSFUN
    tie(g_r[i], g_r[i + threads]) = tagged(r, type);
    tie(g_v[i], g_v[i + threads]) = tagged(v, tag);
#else
    g_r[i] = tagged(r, type);
    g_v[i] = tagged(v, tag);
#endif
    g_image[i] = image;
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 *
 * CUDA execution dimensions must agree with random number generator
 *
 * @param g_v particle velocities (array of size \code{} 2 * nplace \endcode for dsfloat arithmetic)
 * @param g_f particle forces (array of size \code{} nplace \endcode)
 * @param npart number of particles
 * @param nplace number of placeholder particles
 */
template <
    typename vector_type
  , typename vector_type_
  , typename rng_type
  , typename gpu_vector_type
>
__global__ void _finalize(
    float4* g_v
  , gpu_vector_type const* g_f
  , uint npart, uint nplace
  , rng_type rng
)
{
    enum { dimension = vector_type::static_size };

    // read random number generator state from global device memory
    typename rng_type::state_type state = rng[GTID];

    // cache second normal variate for odd dimensions
    bool cached = false;
    typename vector_type::value_type cache;

    // a heat bath collision is performed for all particles of a warp and
    // the collision probability thus requires a correction
    //
    // the probability to have no collision in _all_ threads of a warp
    // must equal the probability to have no collision in a single thread
    // pow(1 - q, warpSize) = 1 - coll_prob_
#ifndef USE_ORIGINAL_ANDERSEN_THERMOSTAT
    float q = 1 - __powf(1 - coll_prob_, 1.f / warpSize);  //< powf() is significantly slower
#endif

    for (uint i = GTID; i < npart; i += GTDIM) {
        // read velocity from global device memory
        unsigned int tag;
        vector_type v;
#ifdef USE_VERLET_DSFUN
        tie(v, tag) = untagged<vector_type>(g_v[i], g_v[i + nplace]);
#else
        tie(v, tag) = untagged<vector_type>(g_v[i]);
#endif

        // is this a deterministic step?
        //
        // to avoid divergent threads within a warp, the stochastic coupling
        // is performed for all particles in the thread as soon as one thread requests it
#ifdef USE_ORIGINAL_ANDERSEN_THERMOSTAT
        if (uniform(rng, state) > coll_prob_) {
#else
        if (__all(uniform(rng, state) > q)) {
#endif
            // read force from global device memory
            vector_type_ f = g_f[i];
            // update velocity
            verlet_kernel::finalize(v, f, timestep_);
        }
        // stochastic coupling with heat bath
        else {
            // parameters for normal distribution
            float const mean = 0;
            float const sigma = sqrt_temperature_;

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

        // write velocity to global device memory
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + nplace]) = tagged(v, tag);
#else
        g_v[i] = tagged(v, tag);
#endif
    }

    // store random number generator state in global device memory
    rng[GTID] = state;
}

} // namespace verlet_nvt_andersen_kernel

template <int dimension, typename rng_type>
verlet_nvt_andersen_wrapper<dimension, rng_type> const
verlet_nvt_andersen_wrapper<dimension, rng_type>::kernel = {
    verlet_nvt_andersen_kernel::timestep_
  , verlet_nvt_andersen_kernel::sqrt_temperature_
  , verlet_nvt_andersen_kernel::coll_prob_
#ifdef USE_VERLET_DSFUN
  , verlet_nvt_andersen_kernel::_integrate<
        fixed_vector<dsfloat, dimension>
      , fixed_vector<float, dimension>
    >
  , verlet_nvt_andersen_kernel::_finalize<
        fixed_vector<dsfloat, dimension>
      , fixed_vector<float, dimension>
      , random::gpu::rand48_rng
    >
#else
  , verlet_nvt_andersen_kernel::_integrate<
        fixed_vector<float, dimension>
      , fixed_vector<float, dimension>
    >
  , verlet_nvt_andersen_kernel::_finalize<
        fixed_vector<float, dimension>
      , fixed_vector<float, dimension>
      , random::gpu::rand48_rng
    >
#endif
};

template class verlet_nvt_andersen_wrapper<3, random::gpu::rand48_rng>;
template class verlet_nvt_andersen_wrapper<2, random::gpu::rand48_rng>;

} // namespace mdsim
} // namespace gpu
} // namespace integrators
} // namespace halmd
