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
#include <numeric>

#include <halmd/mdsim/box.hpp>
#include <halmd/util/log.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim
{

/**
 * Set box edge lengths
 */
template <int dimension, typename float_type>
box<dimension, float_type>::box(options const& vm)
    // dependency injection
    : particle(dynamic_pointer_cast<particle_type>(factory<mdsim::particle<dimension, float_type> >::fetch(vm)))
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
    float_type volume = accumulate(length_.begin(), length_.end(), static_cast<float_type>(1), multiplies<float_type>());
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
    float_type volume = particle->nbox / accumulate(scale_.begin(), scale_.end(), density_, multiplies<float_type>());
    length_ = scale_ * pow(volume, (float_type) 1 / dimension);

    LOG("simulation box edge lengths: " << length_);
    LOG("number density: " << density_);
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class box<3, double>;
template class box<2, double>;
#else
template class box<3, float>;
template class box<2, float>;
#endif

template <int dimension, typename float_type>
boost::shared_ptr<box<dimension, float_type> >
factory<box<dimension, float_type> >::fetch(options const& vm)
{
    if (!box_) {
        box_.reset(new box<dimension, float_type>(vm));
    }
    return box_;
}

#ifndef USE_HOST_SINGLE_PRECISION
template <> factory<box<3, double> >::box_ptr factory<box<3, double> >::box_ = box_ptr();
template class factory<box<3, double> >;
template <> factory<box<2, double> >::box_ptr factory<box<2, double> >::box_ = box_ptr();
template class factory<box<2, double> >;
#else
template <> factory<box<3, float> >::box_ptr factory<box<3, float> >::box_ = box_ptr();
template class factory<box<3, float> >;
template <> factory<box<2, float> >::box_ptr factory<box<2, float> >::box_ = box_ptr();
template class factory<box<2, float> >;
#endif

}} // namespace halmd::mdsim
