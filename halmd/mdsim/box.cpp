/*
 * Copyright © 2010-2011 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <functional>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {

template <int dimension>
box<dimension>::box(matrix_type const& edges)
{
    if (edges.size1() != dimension || edges.size2() != dimension) {
        throw std::invalid_argument("edge vectors have invalid dimensionality");
    }
    edges_ = zero_matrix_type(dimension, dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges_(i, i) = edges(i, i);
    }
    if (norm_inf(edges_ - edges) != 0) {
        throw std::invalid_argument("non-cuboid geomtries are not implemented");
    }
    for (unsigned int i = 0; i < dimension; ++i) {
        length_[i] = edges_(i, i);
    }
    length_half_ = 0.5 * length_;

    LOG("edge lengths of simulation domain: " << length_);
}

template <int dimension>
typename box<dimension>::vector_type
box<dimension>::lowest_corner() const
{
    return -length_half_;
}

template <int dimension>
double box<dimension>::volume() const
{
    return std::accumulate(length_.begin(), length_.end(), 1., std::multiplies<double>());
}

template <typename box_type>
static std::function<typename box_type::vector_type ()>
wrap_lowest_corner(std::shared_ptr<box_type const> self)
{
    return [=]() {
        return self->lowest_corner();
    };
}

template <typename box_type>
static std::function<std::vector<typename box_type::vector_type> ()>
wrap_edges(std::shared_ptr<box_type const> self)
{
    typedef std::vector<typename box_type::vector_type> matrix_type;
    return [=]() -> matrix_type  {
        typename box_type::matrix_type const& input = self->edges();
        matrix_type output(input.size1());
        for (unsigned int i = 0; i < input.size1(); ++i) {
            for (unsigned int j = 0; j < input.size2(); ++j) {
                output[i][j] = input(i, j);
            }
        }
        return output;
    };
}

/**
 * This function is a helper for the H5MD reader, which requires a
 * reference to a value in order to pick the correct data type.
 * Note that the reader does not support matrices, therefore
 * a std::vector of fixed-size vectors is used instead. Both types
 * have the same Lua representation, a table of tables.
 */
template <typename box_type>
static std::function<std::vector<typename box_type::vector_type>& ()>
make_edges()
{
    typedef std::vector<typename box_type::vector_type> matrix_type;
    std::shared_ptr<matrix_type> edges = std::make_shared<matrix_type>();
    return [=]() -> matrix_type& {
        return *edges;
    };
}

template <int dimension>
void box<dimension>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static std::string class_name("box_" + std::to_string(dimension));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<box, std::shared_ptr<box> >(class_name.c_str())
                .def(constructor<matrix_type const&>())
                .property("edges", &wrap_edges<box>)
                .property("lowest_corner", &wrap_lowest_corner<box>)
                .property("length", &box::length)
                .property("volume", &box::volume)
                .scope
                [
                    def("make_edges", &make_edges<box>)
                ]
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
template class box<1>;

} // namespace mdsim
} // namespace halmd
