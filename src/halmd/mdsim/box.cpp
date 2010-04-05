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
#include <cmath>

#include <halmd/mdsim/box.hpp>
#include <halmd/util/log.hpp>

using namespace std;

namespace halmd { namespace mdsim
{

/**
 * Set box edge lengths
 */
template <int dimension, typename float_type>
box<dimension, float_type>::box(particle_ptr const& particle, options const& vm)
    // dependency injection
    : particle(particle)
    // default to cube
    , scale_(1)
{
    // parse options
    if (vm["density"].defaulted() && !vm["box-length"].empty()) {
        length(vm["box-length"].as<float>());
    }
    else {
        density(vm["density"].as<float>());
    }
}

/**
 * Set edge lengths of cuboid
 */
template <int dimension, typename float_type>
void box<dimension, float_type>::length(vector_type const& value)
{
    length_ = value;
    scale_ = length_ / *max_element(length_.begin(), length_.end());
    float_type volume = 1;
    for (size_t i = 0; i < dimension; ++i) {
        volume *= length_[i];
    }
    density_ = particle->nbox / volume;

    LOG("simulation box edge lengths: " << length_);
    LOG("number density: " << density_);
}

/**
 * Set number density
 */
template <int dimension, typename float_type>
void box<dimension, float_type>::density(float_type value)
{
    density_ = value;
    float_type volume = particle->nbox / density_;
    for (size_t i = 0; i < dimension; ++i) {
        volume /= scale_[i];
    }
    length_ = scale_ * pow(volume, (float_type) 1 / dimension);

    LOG("simulation box edge lengths: " << length_);
    LOG("number density: " << density_);
}

// explicit instatiation
#ifndef USE_HOST_SINGLE_PRECISION
template class box<3, double>;
template class box<2, double>;
#else
template class box<3, float>;
template class box<2, float>;
#endif

}} // namespace halmd::mdsim
