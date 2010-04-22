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
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <exception>
#include <numeric>

#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/particle.hpp>
#endif
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/util/logger.hpp>

using namespace boost;
using namespace boost::algorithm;
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
        ntypes.assign(value.begin(), value.end());
    }
    else {
        unsigned int value = vm["particles"].as<unsigned int>();
        if (value < 1) {
            throw logic_error("invalid number of particles");
        }
        ntypes.push_back(value);
    }
    nbox = accumulate(ntypes.begin(), ntypes.end(), 0);
    ntype = ntypes.size();

    vector<string> ntypes_(ntypes.size());
    std::transform(ntypes.begin(), ntypes.end(), ntypes_.begin(), lexical_cast<string, unsigned int>);

    LOG("positional coordinates dimension: " << dimension);
    LOG("number of particles: " << nbox);
    LOG("number of particle types: " << ntype);
    LOG("number of particles per type: " << join(ntypes_, " "));
}

template <int dimension>
typename particle<dimension>::pointer particle<dimension>::create(options const& vm)
{
    LOG_DEBUG("creating PARTICLE");
    if (vm["backend"].as<string>() == "host") {
#ifdef USE_HOST_SINGLE_PRECISION
        return pointer(new host::particle<dimension, float>(vm));
#else
        return pointer(new host::particle<dimension, double>(vm));
#endif
    }
#ifdef WITH_CUDA
    else if (vm["backend"].as<string>() == "gpu_neighbour") {
        return pointer(new gpu::particle<dimension, float>(vm));
    }
#endif
    throw std::runtime_error("not implemented");
}

// explicit instantiation
template class particle<3>;
template class particle<2>;

}} // namespace halmd::mdsim
