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
#include <exception>
#include <numeric>

#include <halmd/algorithm/host/permute.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/util/log.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim { namespace host
{

template <unsigned int dimension, typename float_type>
particle<dimension, float_type>::particle(options const& vm)
    : _Base(vm)
    // allocate particle storage
    , r(nbox)
    , image(nbox)
    , v(nbox)
    , f(nbox)
    , tag(nbox)
    , type(nbox)
    , neighbor(nbox)
{
}

/**
 * Rearrange particles in memory according to an integer index sequence
 *
 * The neighbor lists must be rebuilt after calling this function!
 */
template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::rearrange(std::vector<unsigned int> const& index)
{
    algorithm::host::permute(r.begin(), r.end(), index.begin());
    algorithm::host::permute(image.begin(), image.end(), index.begin());
    algorithm::host::permute(v.begin(), v.end(), index.begin());
    algorithm::host::permute(f.begin(), f.end(), index.begin());
    algorithm::host::permute(tag.begin(), tag.end(), index.begin());
    algorithm::host::permute(type.begin(), type.end(), index.begin());
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class particle<3, double>;
template class particle<2, double>;
#else
template class particle<3, float>;
template class particle<2, float>;
#endif

}}} // namespace halmd::mdsim::host
