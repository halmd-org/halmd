/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <boost/algorithm/string/predicate.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/box.hpp>
#include <halmd/utility/module.hpp>

// workaround
#include <halmd/mdsim/gpu/particle.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu
{

/**
 * Resolve module dependencies
 */
template <int dimension>
void box<dimension>::depends()
{
    modules::required<_Self, particle_type>();
    // workaround
    modules::required<_Self, gpu::particle<dimension, float> >();
}

template <int dimension>
void box<dimension>::select(po::options const& vm)
{
/*    if (!starts_with(vm["backend"].as<string>(), "gpu")) {
        throw unsuitable_module("mismatching option backend");
    }*/
}

/**
 * Set box edge lengths
 */
template <int dimension>
box<dimension>::box(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(modules::fetch<particle_type>(vm))
  // initialize parameters
  , length_half_(.5 * length_)
{}

// explicit instantiation
template class box<3>;
template class box<2>;

}} // namespace mdsim::gpu

template class module<mdsim::gpu::box<3> >;
template class module<mdsim::gpu::box<2> >;

} // namespace halmd
