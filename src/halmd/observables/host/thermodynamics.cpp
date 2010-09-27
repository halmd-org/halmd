/*
 * Copyright © 2010  Felix Höfling
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

#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/utility/module.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace host
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::depends()
{
    modules::depends<_Self, particle_type>::required();
    modules::depends<_Self, force_type>::required();
}

template <int dimension, typename float_type>
thermodynamics<dimension, float_type>::thermodynamics(modules::factory& factory, po::variables_map const& vm)
  : _Base(factory, vm)
  // dependency injection
  , particle(modules::fetch<particle_type>(factory, vm))
  , force(modules::fetch<force_type>(factory, vm))
{
}

template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_kin() const
{
    // compute mean-square velocity
    double vv = 0;
    BOOST_FOREACH(vector_type const& v, particle->v) {
        // assuming unit mass for all particle types
        vv += inner_prod(v, v);
    }
    return .5 * vv / particle->nbox;
}

template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::vector_type thermodynamics<dimension, float_type>::v_cm() const
{
    // compute mean velocity
    vector_type v_cm_(0.);
    BOOST_FOREACH(vector_type const& v, particle->v) {
        v_cm_ += v;
    }
    return v_cm_ / particle->nbox;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class thermodynamics<3, double>;
template class thermodynamics<2, double>;
#else
template class thermodynamics<3, float>;
template class thermodynamics<2, float>;
#endif

}} // namespace observables::host

#ifndef USE_HOST_SINGLE_PRECISION
template class module<observables::host::thermodynamics<3, double> >;
template class module<observables::host::thermodynamics<2, double> >;
#else
template class module<observables::host::thermodynamics<3, float> >;
template class module<observables::host::thermodynamics<2, float> >;
#endif

} // namespace halmd
