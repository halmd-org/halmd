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

#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/host/forces/lj.hpp>
#include <halmd/util/log.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim
{

template <int dimension, typename float_type>
force<dimension, float_type>::force(options const& vm)
    // dependency injection
    : particle(dynamic_pointer_cast<particle_type>(factory<mdsim::particle<dimension, float_type> >::fetch(vm)))
    // allocate result variables
    , virial_(particle->ntype)
{
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class force<3, double>;
template class force<2, double>;
#else
template class force<3, float>;
template class force<2, float>;
#endif

template <int dimension, typename float_type>
boost::shared_ptr<force<dimension, float_type> >
factory<force<dimension, float_type> >::fetch(options const& vm)
{
    if (!force_) {
        force_.reset(new host::forces::lj<dimension, float_type>(vm));
    }
    return force_;
}

#ifndef USE_HOST_SINGLE_PRECISION
template <> factory<force<3, double> >::force_ptr factory<force<3, double> >::force_ = force_ptr();
template class factory<force<3, double> >;
template <> factory<force<2, double> >::force_ptr factory<force<2, double> >::force_ = force_ptr();
template class factory<force<2, double> >;
#else
template <> factory<force<3, float> >::force_ptr factory<force<3, float> >::force_ = force_ptr();
template class factory<force<3, float> >;
template <> factory<force<2, float> >::force_ptr factory<force<2, float> >::force_ = force_ptr();
template class factory<force<2, float> >;
#endif
}} // namespace halmd::mdsim
