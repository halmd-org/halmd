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
        ("box-length,L", po::value<multi_array<float, 1> >(),
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
    register_any_converter<multi_array<float, 1> >();
}

/**
 * Set box edge lengths
 */
template <int dimension>
box<dimension>::box(
    shared_ptr<particle_type> particle
  , vector_type const& length
)
  : length_(length)
  , density_(
        particle->nbox / accumulate(
            length_.begin(), length_.end(), 1., multiplies<double>()
        )
    )
  , length_half_(0.5 * length_)
{
    LOG("edge lengths of simulation box: " << length_);
    LOG("number density: " << density_);
}

/**
 * Set number density
 */
template <int dimension>
box<dimension>::box(
    shared_ptr<particle_type> particle
  , double density
  , vector_type const& ratios
)
  : length_(
        ratios * pow(
            particle->nbox / accumulate(
                ratios.begin(), ratios.end(), density, multiplies<double>()
            )
          , 1. / dimension
        )
    )
  , density_(density)
  , length_half_(0.5 * length_)
{
    LOG("number density: " << density_);
    LOG("edge lengths of simulation box: " << length_);
}

template <typename T>
static void register_lua(lua_State* L, char const* class_name)
{
    typedef typename T::particle_type particle_type;
    typedef typename T::vector_type vector_type;

    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                class_<T, shared_ptr<T> >(class_name)
                    .def(constructor<shared_ptr<particle_type>, vector_type const&>())
                    .def(constructor<shared_ptr<particle_type>, double, vector_type const&>())
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
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        bind(&register_lua<box<3> >, _1, "box_3_")
    ]
    [
        bind(&register_lua<box<2> >, _1, "box_2_")
    ];
}

// explicit instantiation
template class box<3>;
template class box<2>;

} // namespace mdsim

} // namespace halmd
