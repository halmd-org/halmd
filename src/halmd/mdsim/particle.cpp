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

#include <halmd/io/logger.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/particle.hpp>
#endif
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/particle.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Assemble module options
 */
template <int dimension>
void particle<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("backend",
#ifdef WITH_CUDA
         po::value<string>()->default_value("gpu"),
#else
         po::value<string>()->default_value("host"),
#endif
         "computing device type")
        ;

    po::options_description group("Particle");
    group.add_options()
        ("particles,N", po::value<unsigned int>()->default_value(1000),
         "number of particles")
        ("binary,M", po::value<boost::array<unsigned int, 2> >(),
         "binary mixture with A,B particles")
        ;
    desc.add(group);
}

template <int dimension>
particle<dimension>::particle(po::options const& vm)
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

    LOG("MD simulation backend: " << vm["backend"].as<string>());
    LOG("dimension of positional coordinates: " << dimension);
    LOG("number of particles: " << nbox);
    LOG("number of particle types: " << ntype);
    LOG("number of particles per type: " << join(ntypes_, " "));
}

// explicit instantiation
template class particle<3>;
template class particle<2>;

} // namespace mdsim

} // namespace halmd
