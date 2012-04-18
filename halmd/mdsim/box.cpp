/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <cmath>
#include <functional> // std::multiplies
#include <numeric> // std::accumulate

#include <halmd/io/logger.hpp>
#include <halmd/io/utility/hdf5.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {

/**
 * Set box edge lengths
 */
template <int dimension>
box<dimension>::box(vector_type const& length)
  : length_(length)
  , length_half_(0.5 * length_)
{
    LOG("edge lengths of simulation box: " << length_);
}

template <int dimension>
typename box<dimension>::vector_type
box<dimension>::origin() const
{
    return -length_half_;
}

template <int dimension>
vector<typename box<dimension>::vector_type>
box<dimension>::edges() const
{
    vector<vector_type> edges(dimension, 0);
    for (int i = 0; i < dimension; ++i) {
        edges[i][i] = length_[i];
    }
    return edges;
}

template <int dimension>
double box<dimension>::volume() const
{
    return accumulate(length_.begin(), length_.end(), 1., multiplies<double>());
}

template <int dimension>
static int wrap_dimension(box<dimension> const&)
{
    return dimension;
}

template <typename box_type>
static function<typename box_type::vector_type ()>
wrap_origin(shared_ptr<box_type const> box)
{
    return bind(&box_type::origin, box);
}

template <typename box_type>
static function<vector<typename box_type::vector_type> ()>
wrap_edges(boost::shared_ptr<box_type const> box)
{
    return bind(&box_type::edges, box);
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
                .def(constructor<vector_type const&>())
                .property("dimension", &wrap_dimension<dimension>)
                .property("length", &box::length)
                .property("volume", &box::volume)
                .property("origin", &wrap_origin<box>)
                .property("edges", &wrap_edges<box>)
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
