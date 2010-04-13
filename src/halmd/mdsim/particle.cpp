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

#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/util/log.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim
{

template <int dimension>
particle<dimension>::particle(options const& vm)
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
template class particle<3>;
template class particle<2>;

template <int dimension>
typename factory<particle<dimension> >::pointer
factory<particle<dimension> >::fetch(options const& vm)
{
    if (!singleton_) {
#ifdef USE_HOST_SINGLE_PRECISION
        singleton_.reset(new host::particle<dimension, float>(vm));
#else
        singleton_.reset(new host::particle<dimension, double>(vm));
#endif
    }
    return singleton_;
}

template <> factory<particle<3> >::pointer factory<particle<3> >::singleton_ = pointer();
template class factory<particle<3> >;
template <> factory<particle<2> >::pointer factory<particle<2> >::singleton_ = pointer();
template class factory<particle<2> >;

}} // namespace halmd::mdsim
