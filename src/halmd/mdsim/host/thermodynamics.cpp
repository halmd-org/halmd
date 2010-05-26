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

#include <halmd/mdsim/host/thermodynamics.hpp>
#include <halmd/utility/module.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host
{

/**
 * Resolve module dependencies
 */
template <int dimension>
void thermodynamics<dimension>::resolve(po::options const& vm)
{
    module<particle_type>::required(vm);
}

template <int dimension>
thermodynamics<dimension>::thermodynamics(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(module<particle_type>::fetch(vm))
  // allocate result variables
  , virial_(particle->ntype)
{
}

template <int dimension>
double thermodynamics<dimension>::en_kin() const
{
    // compute mean-square velocity
    double vv = 0;
    BOOST_FOREACH(vector_type const& v, particle->v) {
        // assuming unit mass for all particle types
        vv += inner_prod(v, v);
    }
    return .5 * vv / particle->nbox;
}

template <int dimension>
typename thermodynamics<dimension>::vector_type thermodynamics<dimension>::v_cm() const
{
    // compute mean velocity
    vector_type v_cm_(0.);
    BOOST_FOREACH(vector_type const& v, particle->v) {
        v_cm_ += v;
    }
    return v_cm_ / particle->nbox;
}

// explicit instantiation
template class thermodynamics<3>;
template class thermodynamics<2>;

}} // namespace mdsim::host

template class module<mdsim::host::thermodynamics<3> >;
template class module<mdsim::host::thermodynamics<2> >;

} // namespace halmd
