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

#include <cmath>

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
        ("box-length,L", po::value<multi_array<float, 1> >()->conflicts("density"),
         "edge lengths of simulation box")
        ("box-ratios", po::value<multi_array<float, 1> >()->conflicts("box-length"),
         "aspect ratios of simulation box (specify relative edge lengths)")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
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
  , length_half_(0.5 * length_)
{
    density_ = particle->nbox / volume();

    LOG("number density: " << density_);
    LOG("edge lengths of simulation box: " << length_);
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
  : density_(density)
{
    double volume = particle->nbox / density;
    double det = accumulate(ratios.begin(), ratios.end(), 1., multiplies<double>());
    length_ = ratios * pow(volume / det, 1. / dimension);
    length_half_ = .5 * length_;

    LOG("number density: " << density_);
    LOG("edge lengths of simulation box: " << length_);
}

template <int dimension>
void box<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    string class_name("box_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                class_<box, shared_ptr<box> >(class_name.c_str())
                    .def(constructor<shared_ptr<particle_type>, vector_type const&>())
                    .def(constructor<shared_ptr<particle_type>, double, vector_type const&>())
                    .scope
                    [
                        def("options", &box::options)
                    ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &box<3>::luaopen
    ]
    [
        &box<2>::luaopen
    ];
}

// explicit instantiation
template class box<3>;
template class box<2>;

} // namespace mdsim

} // namespace halmd
