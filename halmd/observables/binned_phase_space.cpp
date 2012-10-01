/*
 * Copyright © 2011  Felix Höfling
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

#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

#include <halmd/observables/binned_phase_space.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {

template <int dimension>
typename binned_phase_space<dimension>::position_grid_type
binned_phase_space<dimension>::position() const
{
    // store references on position grid per axis
    array<vector<double> const*, dimension> position_ref;
    for (unsigned int i = 0; i < dimension; ++i) {
        position_ref[i] = &position(i);
    }

    // determine grid extents
    array<size_t, dimension + 1> nbin;
    for (unsigned int i = 0; i < dimension; ++i) {
        nbin[i] = position_ref[i]->size();
    }
    nbin[dimension] = dimension; // the last dimension refers to the coordinates of a position vector

    // allocate memory
    position_grid_type grid(nbin);

    // iterate over all elements of grid
    assert(grid.num_elements() % dimension == 0);
    double *end = grid.data() + grid.num_elements();
    for (double *it = grid.data(); it != end; ) {
        ssize_t n = (it - grid.data()) / dimension; //< 1-dimensional, global index in cell grid
        // construct d-dimensional grid index
        array<size_t, dimension> idx;
        for (int i = dimension - 1; i >= 0; --i) {
            idx[i] = n % nbin[i];
            n /= nbin[i];
        }
        // convert to position coordinates of grid point
        for (int i = 0; i < dimension; ++i) {
            *it++ = position_ref[i]->at(idx[i]);
        }
    }

    return grid;
}

template <typename binned_phase_space_type>
static typename binned_phase_space_type::slot_function_type
acquire_wrapper(shared_ptr<binned_phase_space_type> self)
{
    return bind(&binned_phase_space_type::acquire, self);
}

template <typename binned_phase_space_type>
static function<vector<double> const& ()>
wrap_position(shared_ptr<binned_phase_space_type> self, unsigned int axis)
{
    return bind(&binned_phase_space_type::position, self, axis);
}

template <typename binned_phase_space_type>
static function<typename binned_phase_space_type::position_grid_type ()>
wrap_position_grid(shared_ptr<binned_phase_space_type> self)
{
    return bind(&binned_phase_space_type::position, self);
}

template <int dimension>
void binned_phase_space<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("binned_phase_space_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<binned_phase_space, shared_ptr<binned_phase_space> >(class_name.c_str())
                .property("acquire", &acquire_wrapper<binned_phase_space>)
                .def("position", &wrap_position_grid<binned_phase_space>)
                .def("position", &wrap_position<binned_phase_space>)
                .def("on_acquire", &binned_phase_space::on_acquire)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_binned_phase_space(lua_State* L)
{
    binned_phase_space<2>::luaopen(L);
    binned_phase_space<3>::luaopen(L);
    return 0;
}

// explicit instantiation
template class binned_phase_space<2>;
template class binned_phase_space<3>;

} // namespace observables
} // namespace halmd
