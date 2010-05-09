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

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/velocity/file.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace velocity
{

using namespace boost;
using namespace std;

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void file<dimension, float_type>::resolve(po::options const& vm)
{
    if (!vm.count("trajectory-sample")) {
        throw module_error("inept module " + module<file>::name());
    }

    module<reader_type>::required(vm);
    module<sample_type>::required(vm);
    module<particle_type>::required(vm);
}

template <int dimension, typename float_type>
file<dimension, float_type>::file(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , reader(module<reader_type>::fetch(vm))
  , sample(module<sample_type>::fetch(vm))
  , particle(module<particle_type>::fetch(vm))
{}

/**
 * set particle velocities
 */
template <int dimension, typename float_type>
void file<dimension, float_type>::set()
{
    for (size_t j = 0, i = 0; j < particle->ntype; i += particle->ntypes[j], ++j) {
        copy(sample->v[j]->begin(), sample->v[j]->end(), &particle->v[i]);
    }

    LOG("set particle velocities from trajectory sample");
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class file<3, double>;
template class file<2, double>;
#else
template class file<3, float>;
template class file<2, float>;
#endif

}}} // namespace mdsim::host::velocity

#ifndef USE_HOST_SINGLE_PRECISION
template class module<mdsim::host::velocity::file<3, double> >;
template class module<mdsim::host::velocity::file<2, double> >;
#else
template class module<mdsim::host::velocity::file<3, float> >;
template class module<mdsim::host::velocity::file<2, float> >;
#endif

} // namespace halmd
