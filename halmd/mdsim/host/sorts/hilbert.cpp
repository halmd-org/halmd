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
#include <cmath>

#include <halmd/mdsim/host/sorts/hilbert.hpp>
#include <halmd/mdsim/sorts/hilbert_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace sorts {

template <int dimension, typename float_type>
hilbert<dimension, float_type>::hilbert(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<binning_type> binning
  , std::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , binning_(binning)
  , logger_(logger)
{
    cell_size_type const& ncell = binning_->ncell();
    vector_type const& cell_length = binning_->cell_length();
    cache_proxy<cell_array_type const> cell = binning_->cell();

    // set Hilbert space-filling curve recursion depth
    unsigned int ncell_max = *max_element(ncell.begin(), ncell.end());
    unsigned int depth = static_cast<int>(std::ceil(std::log(static_cast<double>(ncell_max)) / M_LN2));
    // 32-bit integer for 2D Hilbert code allows a maximum of 16/10 levels
    depth = min((dimension == 3) ? 10U : 16U, depth);

    LOG("vertex recursion depth: " << depth);

    // generate 1-dimensional Hilbert curve mapping of cell lists
    typedef std::pair<cell_list const*, unsigned int> pair;
    std::vector<pair> pairs;
    cell_size_type x;
    for (x[0] = 0; x[0] < ncell[0]; ++x[0]) {
        for (x[1] = 0; x[1] < ncell[1]; ++x[1]) {
            if (dimension == 3) {
                for (x[2] = 0; x[2] < ncell[2]; ++x[2]) {
                    vector_type r(x);
                    r = element_prod(r + vector_type(0.5), cell_length);
                    pairs.push_back(std::make_pair(&(*cell)(x), map(r, depth)));
                }
            }
            else {
                vector_type r(x);
                r = element_prod(r + vector_type(0.5), cell_length);
                pairs.push_back(std::make_pair(&(*cell)(x), map(r, depth)));
            }
        }
    }
    stable_sort(pairs.begin(), pairs.end(), bind(&pair::second, _1) < bind(&pair::second, _2));
    map_.clear();
    map_.reserve(cell->size());
    transform(pairs.begin(), pairs.end(), back_inserter(map_), bind(&pair::first, _1));
}

/**
 * Order particles after Hilbert space-filling curve
 */
template <int dimension, typename float_type>
void hilbert<dimension, float_type>::order()
{
    LOG_TRACE("order particles after Hilbert space-filling curve");
    {
        scoped_timer_type timer(runtime_.order);
        std::vector<unsigned int> index;
        {
            scoped_timer_type timer(runtime_.map);
            // particle binning
            binning_->cell();
            // generate index sequence according to Hilbert-sorted cells
            index.reserve(particle_->nparticle());
            for (cell_list const* cell : map_) {
                for (unsigned int p : *cell) {
                    index.push_back(p);
                }
            }
        }

        // reorder particles in memory
        particle_->rearrange(index);
    }
    on_order_();
}

/**
 * Map 3-/2-dimensional point to 1-dimensional point on Hilbert space curve
 */
template <int dimension, typename float_type>
unsigned int hilbert<dimension, float_type>::map(vector_type r, unsigned int depth)
{
    r = element_div(r, static_cast<vector_type>(box_->length()));

    return mdsim::sorts::hilbert_kernel::map(r, depth);
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(hilbert<dimension, float_type> const&)
{
    return hilbert<dimension, float_type>::module_name();
}

template <typename sort_type>
static std::function<void ()>
wrap_order(std::shared_ptr<sort_type> self)
{
    return [=]() {
        self->order();
    };
}

template <int dimension, typename float_type>
void hilbert<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static string class_name("hilbert_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("sorts")
                [
                    class_<hilbert, std::shared_ptr<hilbert> >(class_name.c_str())
                        .def(constructor<
                            std::shared_ptr<particle_type>
                          , std::shared_ptr<box_type const>
                          , std::shared_ptr<binning_type>
                          , std::shared_ptr<logger_type>
                        >())
                        .property("module_name", &module_name_wrapper<dimension, float_type>)
                        .property("order", &wrap_order<hilbert>)
                        .def("on_order", &hilbert::on_order)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("order", &runtime::order)
                                .def_readonly("map", &runtime::map)
                        ]
                        .def_readonly("runtime", &hilbert::runtime_)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_sorts_hilbert(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    hilbert<3, double>::luaopen(L);
    hilbert<2, double>::luaopen(L);
#else
    hilbert<3, float>::luaopen(L);
    hilbert<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class hilbert<3, double>;
template class hilbert<2, double>;
#else
template class hilbert<3, float>;
template class hilbert<2, float>;
#endif

} // namespace sorts
} // namespace host
} // namespace mdsim
} // namespace halmd
