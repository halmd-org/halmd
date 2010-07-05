/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/velocity/boltzmann.hpp>
#include <halmd/mdsim/gpu/velocity/boltzmann_kernel.cuh>
#include <halmd/utility/module.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace velocity
{

using namespace boost;
using namespace std;

/**
 * Assemble module options
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void boltzmann<dimension, float_type, RandomNumberGenerator>::options(po::options_description& desc)
{
    po::options_description group("Boltzmann velocity distribution");
    group.add_options()
        ("temperature,K", po::value<float>()->default_value(1.12),
         "Boltzmann distribution temperature")
        ;
    desc.add(group);
}

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void boltzmann<dimension, float_type, RandomNumberGenerator>::depends()
{
    modules::depends<_Self, particle_type>::required();
    modules::depends<_Self, random_type>::required();
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void boltzmann<dimension, float_type, RandomNumberGenerator>::select(po::options const& vm)
{
    if (vm["velocity"].as<string>() != "boltzmann") {
        throw unsuitable_module("mismatching option velocity");
    }
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
boltzmann<dimension, float_type, RandomNumberGenerator>::boltzmann(modules::factory& factory, po::options const& vm)
  : _Base(factory, vm)
  // dependency injection
  , particle(modules::fetch<particle_type>(factory, vm))
  , random(modules::fetch<random_type>(factory, vm))
{
    // parse options
    temp_ = vm["temperature"].as<float>();
}

/**
 * Initialise velocities from Maxwell-Boltzmann distribution
 *
 * The particle velocities need to fullfill two constraints:
 *
 *  1. center of mass velocity shall be zero
 *  2. temperature of the distribution shall equal exactly the given value
 *
 * The above order is chosen as shifting the center of mass velocity
 * means altering the first moment of the velocity distribution, which
 * in consequence affects the second moment, i.e. the temperature.
 *
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void boltzmann<dimension, float_type, RandomNumberGenerator>::set()
{
    LOG("assigning Maxwell-Boltzmann velocity distribution: T = " << temp_);

#ifdef USE_VERLET_DSFUN
    cuda::memset(particle->g_v, 0, particle->g_v.capacity());
#endif

    typedef boltzmann_wrapper<dimension> _gpu;

    /** block sum of velocity */
    cuda::vector<gpu_vector_type> g_vcm;
    /** block sum of squared velocity */
    cuda::vector<dsfloat> g_vv;
#ifdef USE_VERLET_DSFUN
    g_vcm.resize(2 * _gpu::BLOCKS);
#else
    g_vcm.resize(_gpu::BLOCKS);
#endif
    g_vv.resize(_gpu::BLOCKS);

    // generate Maxwell-Boltzmann distributed velocities,
    // assuming equal (unit) mass for all particle types
    // random->normal<dimension>(particle->g_v, temp_);

    // determine g_vcm and g_vv

    // set center of mass velocity to zero and reduce squared velocity
//     cuda::configure(_gpu::BLOCKS, _gpu::THREADS);
    cuda::configure(particle->dim.grid, particle->dim.block);
    _gpu::shift_velocity(particle->g_v, particle->nbox, particle->dim.threads(), g_vcm, g_vv);

    // rescale velocities to accurate temperature
//     cuda::configure(_gpu::BLOCKS, _gpu::THREADS);
    cuda::configure(particle->dim.grid, particle->dim.block);
    _gpu::scale_velocity(particle->g_v, particle->nbox, particle->dim.threads(), g_vv, temp_);

/*
    // center velocities around origin, then rescale to exactly
    // match the desired temperature;
    // temp = vv / dimension
    // vv changes to vv - v_cm^2 after shifting
    float_type scale = sqrt(temp_ * dimension / (vv - inner_prod(v_cm, v_cm)));
    shift_rescale(-v_cm, scale);*/


#ifdef USE_HILBERT_ORDER
    // make thermostat independent of neighbour list update frequency or skin
//     order_velocities(); boltzmann is not a thermostat!
#endif



//    LOG_DEBUG("velocities rescaled by factor " << scale);
}

/**
 * Assign new velocities from Gaussian distribution
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
pair<typename boltzmann<dimension, float_type, RandomNumberGenerator>::vector_type, float_type>
inline boltzmann<dimension, float_type, RandomNumberGenerator>::gaussian(float_type sigma)
{
    vector_type v_cm = 0;
    float_type vv = 0;
 /*   float_type r;
    bool r_valid = false;

    BOOST_FOREACH (vector_type& v, particle->v) {
        // assign two components at a time
        for (unsigned i=0; i < dimension-1; i+=2) {
            random->normal(v[i], v[i+1], sigma);
        }
        // handle last component separately for odd dimensions
        if (dimension % 2 == 1) {
            if (r_valid) {
                v[dimension-1] = r;
            }
            else {
                random->normal(v[dimension-1], r, sigma);
            }
            r_valid = !r_valid;
        }
        v_cm += v;
        vv += inner_prod(v, v);
    }

    v_cm /= particle->v.size();
    vv /= particle->v.size();*/
    return make_pair(v_cm, vv);
}

/**
 * Shift all velocities by 'v'
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
inline void boltzmann<dimension, float_type, RandomNumberGenerator>::shift(vector_type const& v_shift)
{
/*    BOOST_FOREACH (vector_type& v, particle->v) {
        v += v_shift;
    }*/
}

/**
 * Rescale magnitude of all velocities by factor 'scale'
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
inline void boltzmann<dimension, float_type, RandomNumberGenerator>::rescale(float_type scale)
{
//     BOOST_FOREACH (vector_type& v, particle->v) {
//         v *= scale;
//     }
//     LOG("velocities rescaled by factor " << scale);
}

/**
 * First shift, then rescale all velocities
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
inline void boltzmann<dimension, float_type, RandomNumberGenerator>::shift_rescale(
    vector_type const& v_shift,
    float_type scale)
{
//     BOOST_FOREACH (vector_type& v, particle->v) {
//         v += v_shift;
//         v *= scale;
//     }
}

}}} // namespace mdsim::gpu::velocity

// explicit instantiation
template class mdsim::gpu::velocity::boltzmann<3, float, random::gpu::rand48>;
template class mdsim::gpu::velocity::boltzmann<2, float, random::gpu::rand48>;

template class module<mdsim::gpu::velocity::boltzmann<3, float, random::gpu::rand48> >;
template class module<mdsim::gpu::velocity::boltzmann<2, float, random::gpu::rand48> >;

} // namespace halmd
