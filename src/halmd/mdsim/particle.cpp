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
#include <boost/array.hpp>
#include <exception>
#include <numeric>

#include <halmd/mdsim/particle.hpp>
#include <halmd/util/log.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim
{

template <int dimension, typename float_type>
particle<dimension, float_type>::particle(options const& vm)
{
    // parse options
    if (!vm["binary"].empty() && vm["particles"].defaulted()) {
        array<unsigned int, 2> value = vm["binary"].as<array<unsigned int, 2> >();
        if (*min_element(value.begin(), value.end()) < 1) {
            throw logic_error("invalid number of A or B particles");
        }
        nbox = accumulate(value.begin(), value.end(), 0);
        ntype = 2;
    }
    else {
        unsigned int value = vm["particles"].as<unsigned int>();
        if (value < 1) {
            throw logic_error("invalid number of particles");
        }
        nbox = value;
        ntype = 1;
    }

    LOG("positional coordinates dimension: " << dimension);
    LOG("number of particles: " << nbox);
    LOG("number of particle types: " << ntype);
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class particle<3, double>;
template class particle<2, double>;
#else
template class particle<3, float>;
template class particle<2, float>;
#endif

}} // namespace halmd::mdsim
