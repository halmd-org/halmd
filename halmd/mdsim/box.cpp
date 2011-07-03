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
#include <halmd/io/utility/hdf5.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Set box edge lengths
 */
template <int dimension>
box<dimension>::box(
    size_t nbox
  , vector_type const& length
)
  : length_(length)
  , length_half_(0.5 * length_)
{
    density_ = nbox / volume();

    LOG("total number of particles: " << nbox);
    LOG("number density: " << density_);
    LOG("edge lengths of simulation box: " << length_);
}

/**
 * Set number density
 */
template <int dimension>
box<dimension>::box(
    size_t nbox
  , double density
  , vector_type const& ratios
)
  : density_(density)
{
    double volume = nbox / density;
    double det = accumulate(ratios.begin(), ratios.end(), 1., multiplies<double>());
    length_ = ratios * pow(volume / det, 1. / dimension);
    length_half_ = .5 * length_;

    LOG("total number of particles: " << nbox);
    LOG("number density: " << density_);
    LOG("edge lengths of simulation box: " << length_);
}

template <int dimension>
static int wrap_dimension(box<dimension> const&)
{
    return dimension;
}

template <int dimension>
void box<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("box_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<box, shared_ptr<box> >(class_name.c_str())
                .def(constructor<size_t, vector_type const&>())
                .def(constructor<size_t, double, vector_type const&>())
                .property("dimension", &wrap_dimension<dimension>)
                .property("length", &box::length)
                .property("density", &box::density)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_box(lua_State* L)
{
    box<3>::luaopen(L);
    box<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class box<3>;
template class box<2>;

} // namespace mdsim

} // namespace halmd
