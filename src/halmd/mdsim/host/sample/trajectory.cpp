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

#include <halmd/mdsim/host/sample/trajectory.hpp>
#include <halmd/util/logger.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace sample
{

template <int dimension, typename float_type>
trajectory<dimension, float_type>::trajectory(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(module<particle_type>::fetch(vm))
  , box(module<box_type>::fetch(vm))
{}

/**
 * Sample trajectory
 */
template <int dimension, typename float_type>
void trajectory<dimension, float_type>::acquire()
{
    r.resize(particle->ntype);
    v.resize(particle->ntype);
    for (size_t i = 0; i < particle->ntype; ++i) {
        r[i].reset(new position_sample_vector(particle->ntypes[i]));
        v[i].reset(new velocity_sample_vector(particle->ntypes[i]));
    }
    for (size_t i = 0; i < particle->nbox; ++i) {
        // periodically extended particle position
        (*r[particle->type[i]])[particle->tag[i]] = particle->r[i] + element_prod(particle->image[i], box->length());
        // particle velocity
        (*v[particle->type[i]])[particle->tag[i]] = particle->v[i];
    }
}

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void trajectory<dimension, float_type>::resolve(po::options const& vm)
{
    module<particle_type>::resolve(vm);
    module<box_type>::resolve(vm);
}

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
po::options_description
trajectory<dimension, float_type>::options()
{
    po::options_description desc;
    desc.add_options()
        // no module options
        ;
    return desc;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class trajectory<3, double>;
template class trajectory<2, double>;
#else
template class trajectory<3, float>;
template class trajectory<2, float>;
#endif

}}} // namespace mdsim::host::sample

#ifndef USE_HOST_SINGLE_PRECISION
template class module<mdsim::host::sample::trajectory<3, double> >;
template class module<mdsim::host::sample::trajectory<2, double> >;
#else
template class module<mdsim::host::sample::trajectory<3, float> >;
template class module<mdsim::host::sample::trajectory<2, float> >;
#endif

} // namespace halmd

