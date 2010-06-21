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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/velocity/boltzmann.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace velocity
{

using namespace boost;
using namespace std;

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
void boltzmann<dimension, float_type>::options(po::options_description& desc)
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
template <int dimension, typename float_type>
void boltzmann<dimension, float_type>::depends()
{
    modules::depends<_Self, particle_type>::required();
    modules::depends<_Self, random_type>::required();
}

template <int dimension, typename float_type>
void boltzmann<dimension, float_type>::select(po::options const& vm)
{
    if (vm["velocity"].as<string>() != "boltzmann") {
        throw unsuitable_module("mismatching option velocity");
    }
}

template <int dimension, typename float_type>
boltzmann<dimension, float_type>::boltzmann(modules::factory& factory, po::options const& vm)
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
 */
template <int dimension, typename float_type>
void boltzmann<dimension, float_type>::set()
{
    pair<vector_type, float_type> p;
    vector_type& v_cm = p.first;
    float_type& vv = p.second;

    // assuming equal (unit) mass for all particle types
    p = gaussian(sqrt(temp_));

    // center velocities around origin, then rescale to exactly
    // match the desired temperature;
    // temp = vv / dimension
    // vv changes to vv - v_cm^2 after shifting
    float_type scale = sqrt(temp_ * dimension / (vv - inner_prod(v_cm, v_cm)));
    shift_rescale(-v_cm, scale);

    LOG_DEBUG("velocities rescaled by factor " << scale);
    LOG("assigned Maxwell-Boltzmann velocity distribution: T = " << temp_);
}

/**
 * Assign new velocities from Gaussian distribution
 */
template <int dimension, typename float_type>
pair<typename boltzmann<dimension, float_type>::vector_type, float_type>
inline boltzmann<dimension, float_type>::gaussian(float_type sigma)
{
    vector_type v_cm = 0;
    float_type vv = 0;
    float_type r;
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
    vv /= particle->v.size();
    return make_pair(v_cm, vv);
}

/**
 * Shift all velocities by 'v'
 */
template <int dimension, typename float_type>
inline void boltzmann<dimension, float_type>::shift(vector_type const& v_shift)
{
    BOOST_FOREACH (vector_type& v, particle->v) {
        v += v_shift;
    }
}

/**
 * Rescale magnitude of all velocities by factor 'scale'
 */
template <int dimension, typename float_type>
inline void boltzmann<dimension, float_type>::rescale(float_type scale)
{
    BOOST_FOREACH (vector_type& v, particle->v) {
        v *= scale;
    }
    LOG("velocities rescaled by factor " << scale);
}

/**
 * First shift, then rescale all velocities
 */
template <int dimension, typename float_type>
inline void boltzmann<dimension, float_type>::shift_rescale(
    vector_type const& v_shift,
    float_type scale)
{
    BOOST_FOREACH (vector_type& v, particle->v) {
        v += v_shift;
        v *= scale;
    }
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class boltzmann<3, double>;
template class boltzmann<2, double>;
#else
template class boltzmann<3, float>;
template class boltzmann<2, float>;
#endif

}}} // namespace mdsim::host::velocity

#ifndef USE_HOST_SINGLE_PRECISION
template class module<mdsim::host::velocity::boltzmann<3, double> >;
template class module<mdsim::host::velocity::boltzmann<2, double> >;
#else
template class module<mdsim::host::velocity::boltzmann<3, float> >;
template class module<mdsim::host::velocity::boltzmann<2, float> >;
#endif

} // namespace halmd
