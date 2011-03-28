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
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <cmath>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/sorts/hilbert.hpp>
#include <halmd/mdsim/sorts/hilbert_kernel.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace sorts
{

template <int dimension, typename float_type>
hilbert<dimension, float_type>::hilbert(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<neighbour_type> neighbour
)
  // dependency injection
  : particle(particle)
  , box(box)
  , neighbour(neighbour)
{
    // set Hilbert space-filling curve recursion depth
    unsigned int ncell = *max_element(neighbour->ncell_.begin(), neighbour->ncell_.end());
    unsigned int depth = static_cast<int>(std::ceil(std::log(static_cast<double>(ncell)) / M_LN2));
    // 32-bit integer for 2D Hilbert code allows a maximum of 16/10 levels
    depth = min((dimension == 3) ? 10U : 16U, depth);

    LOG("Hilbert vertex recursion depth: " << depth);

    // generate 1-dimensional Hilbert curve mapping of cell lists
    typedef std::pair<cell_list*, unsigned int> pair;
    std::vector<pair> pairs;
    cell_size_type x;
    for (x[0] = 0; x[0] < neighbour->ncell_[0]; ++x[0]) {
        for (x[1] = 0; x[1] < neighbour->ncell_[1]; ++x[1]) {
            if (dimension == 3) {
                for (x[2] = 0; x[2] < neighbour->ncell_[2]; ++x[2]) {
                    vector_type r(x);
                    r = element_prod(r + vector_type(0.5), neighbour->cell_length_);
                    pairs.push_back(std::make_pair(&neighbour->cell_(x), map(r, depth)));
                }
            }
            else {
                vector_type r(x);
                r = element_prod(r + vector_type(0.5), neighbour->cell_length_);
                pairs.push_back(std::make_pair(&neighbour->cell_(x), map(r, depth)));
            }
        }
    }
    stable_sort(pairs.begin(), pairs.end(), bind(&pair::second, _1) < bind(&pair::second, _2));
    cell_.clear();
    cell_.reserve(neighbour->cell_.size());
    transform(pairs.begin(), pairs.end(), back_inserter(cell_), bind(&pair::first, _1));
}

/**
 * Order particles after Hilbert space-filling curve
 */
template <int dimension, typename float_type>
void hilbert<dimension, float_type>::order()
{
    // particle binning
    neighbour->update_cells();
    // generate index sequence according to Hilbert-sorted cells
    std::vector<unsigned int> index;
    index.reserve(particle->nbox);
    BOOST_FOREACH(cell_list* cell, cell_) {
        BOOST_FOREACH(unsigned int p, *cell) {
            index.push_back(p);
        }
    }
    // reorder particles in memory
    particle->rearrange(index);
}

/**
 * Map 3-/2-dimensional point to 1-dimensional point on Hilbert space curve
 */
template <int dimension, typename float_type>
unsigned int hilbert<dimension, float_type>::map(vector_type r, unsigned int depth)
{
    r = element_div(r, static_cast<vector_type>(box->length()));

    return mdsim::sorts::hilbert_kernel::map(r, depth);
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(hilbert<dimension, float_type> const&)
{
    return hilbert<dimension, float_type>::module_name();
}

template <int dimension, typename float_type>
void hilbert<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("hilbert_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    namespace_("sorts")
                    [
                        class_<hilbert, shared_ptr<_Base>, _Base>(class_name.c_str())
                            .def(constructor<
                                shared_ptr<particle_type>
                              , shared_ptr<box_type>
                              , shared_ptr<neighbour_type>
                            >())
                            .property("module_name", &module_name_wrapper<dimension, float_type>)
                    ]
                ]
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &hilbert<3, double>::luaopen
    ]
    [
        &hilbert<2, double>::luaopen
    ];
#else
    [
        &hilbert<3, float>::luaopen
    ]
    [
        &hilbert<2, float>::luaopen
    ];
#endif
}

} // namespace

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class hilbert<3, double>;
template class hilbert<2, double>;
#else
template class hilbert<3, float>;
template class hilbert<2, float>;
#endif

}}} // namespace mdsim::host::sorts

} // namespace halmd
