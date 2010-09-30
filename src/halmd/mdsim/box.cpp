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

#include <algorithm>
#include <cmath>
#include <numeric>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Assemble module options
 */
template <int dimension>
void box<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("density,d", po::value<float>()->default_value(0.75),
         "particle density")
        ("box-length,L", po::value<float>(),
         "simulation box length")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    using namespace lua_wrapper;
    register_any_converter<float>();
}


/**
 * Resolve module dependencies
 */
template <int dimension>
void box<dimension>::depends()
{
    modules::depends<_Self, particle_type>::required();
}

/**
 * Set box edge lengths
 */
template <int dimension>
box<dimension>::box(modules::factory& factory, po::variables_map const& vm)
  // dependency injection
  : particle(modules::fetch<particle_type>(factory, vm))
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
template <int dimension>
void box<dimension>::length(vector_type const& value)
{
    length_ = value;
    length_half_ = .5 * length_;
    scale_ = length_ / *max_element(length_.begin(), length_.end());
    double volume = accumulate(length_.begin(), length_.end(), 1., multiplies<double>());
    density_ = particle->nbox / volume;

    LOG("edge lengths of simulation box: " << length_);
    LOG("number density: " << density_);
}

/**
 * Set number density
 */
template <int dimension>
void box<dimension>::density(double value)
{
    density_ = value;
    double volume = particle->nbox / accumulate(scale_.begin(), scale_.end(), density_, multiplies<double>());
    length_ = scale_ * pow(volume, 1. / dimension);
    length_half_ = .5 * length_;

    LOG("simulation box edge lengths: " << length_);
    LOG("number density: " << density_);
}

template <typename T>
static void register_lua(char const* class_name)
{
    using namespace luabind;
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                class_<T, shared_ptr<T> >(class_name)
                    .scope
                    [
                        def("options", &T::options)
                    ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    register_lua<box<3> >("box_3_");
    register_lua<box<2> >("box_2_");
}

// explicit instantiation
template class box<3>;
template class box<2>;

} // namespace mdsim

template class module<mdsim::box<3> >;
template class module<mdsim::box<2> >;

} // namespace halmd
