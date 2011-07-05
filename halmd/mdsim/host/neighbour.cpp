/*
 * Copyright © 2008-2011  Peter Colberg
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

#include <boost/foreach.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/neighbour.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {

/**
 * construct neighbour list module
 *
 * @param particle mdsim::host::particle instance
 * @param box mdsim::box instance
 * @param cutoff force cutoff radius
 * @param skin neighbour list skin
 */
template <int dimension, typename float_type>
neighbour<dimension, float_type>::neighbour(
    shared_ptr<particle_type const> particle
  , shared_ptr<box_type const> box
  , shared_ptr<binning_type const> binning
  , matrix_type const& r_cut
  , double skin
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , binning_(binning)
  // allocate parameters
  , neighbour_(particle_->nbox)
  , r_skin_(skin)
  , rr_cut_skin_(particle_->ntype, particle_->ntype)
{
    matrix_type r_cut_skin(particle_->ntype, particle_->ntype);
    typename matrix_type::value_type r_cut_max = 0;
    for (size_t i = 0; i < particle_->ntype; ++i) {
        for (size_t j = i; j < particle_->ntype; ++j) {
            r_cut_skin(i, j) = r_cut(i, j) + r_skin_;
            rr_cut_skin_(i, j) = std::pow(r_cut_skin(i, j), 2);
            r_cut_max = max(r_cut_skin(i, j), r_cut_max);
        }
    }

    LOG("neighbour list skin: " << r_skin_);
}

/**
 * Update neighbour lists
 */
template <int dimension, typename float_type>
void neighbour<dimension, float_type>::update()
{
    cell_size_type const& ncell = binning_->ncell();
    cell_size_type i;
    for (i[0] = 0; i[0] < ncell[0]; ++i[0]) {
        for (i[1] = 0; i[1] < ncell[1]; ++i[1]) {
            if (dimension == 3) {
                for (i[2] = 0; i[2] < ncell[2]; ++i[2]) {
                    update_cell_neighbours(i);
                }
            }
            else {
                update_cell_neighbours(i);
            }
        }
    }
}

/**
 * Update neighbour lists for a single cell
 */
template <int dimension, typename float_type>
void neighbour<dimension, float_type>::update_cell_neighbours(cell_size_type const& i)
{
    cell_lists const& cell = binning_->cell();
    cell_size_type const& ncell = binning_->ncell();

    BOOST_FOREACH(size_t p, cell(i)) {
        // empty neighbour list of particle
        neighbour_[p].clear();

        cell_diff_type j;
        for (j[0] = -1; j[0] <= 1; ++j[0]) {
            for (j[1] = -1; j[1] <= 1; ++j[1]) {
                if (dimension == 3) {
                    for (j[2] = -1; j[2] <= 1; ++j[2]) {
                        // visit half of 26 neighbour cells due to pair potential
                        if (j[0] == 0 && j[1] == 0 && j[2] == 0) {
                            goto self;
                        }
                        // update neighbour list of particle
                        cell_size_type k = element_mod(static_cast<cell_size_type>(static_cast<cell_diff_type>(i + ncell) + j), ncell);
                        compute_cell_neighbours<false>(p, cell(k));
                    }
                }
                else {
                    // visit half of 8 neighbour cells due to pair potential
                    if (j[0] == 0 && j[1] == 0) {
                        goto self;
                    }
                    // update neighbour list of particle
                    cell_size_type k = element_mod(static_cast<cell_size_type>(static_cast<cell_diff_type>(i + ncell) + j), ncell);
                    compute_cell_neighbours<false>(p, cell(k));
                }
            }
        }
self:
        // visit this cell
        compute_cell_neighbours<true>(p, cell(i));
    }
}

/**
 * Update neighbour list of particle
 */
template <int dimension, typename float_type>
template <bool same_cell>
void neighbour<dimension, float_type>::compute_cell_neighbours(size_t i, cell_list const& c)
{
    BOOST_FOREACH(size_t j, c) {
        // skip identical particle and particle pair permutations if same cell
        if (same_cell
         && particle_->type[j] <= particle_->type[i] //< lexical order of (type, tag)
         && particle_->tag[j] <= particle_->tag[i]
        ) {
            continue;
        }

        // particle distance vector
        vector_type r = particle_->r[i] - particle_->r[j];
        box_->reduce_periodic(r);
        // particle types
        size_t a = particle_->type[i];
        size_t b = particle_->type[j];
        // squared particle distance
        float_type rr = inner_prod(r, r);

        // enforce cutoff radius with neighbour list skin
        if (rr >= rr_cut_skin_(a, b)) {
            continue;
        }

        // add particle to neighbour list
        neighbour_[i].push_back(j);
    }
}

template <int dimension, typename float_type>
void neighbour<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("neighbour_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                class_<neighbour, shared_ptr<mdsim::neighbour>, mdsim::neighbour>(class_name.c_str())
                    .def(constructor<
                         shared_ptr<particle_type const>
                       , shared_ptr<box_type const>
                       , shared_ptr<binning_type const>
                       , matrix_type const&
                       , double
                     >())
                    .property("r_skin", &neighbour::r_skin)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_neighbour(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    neighbour<3, double>::luaopen(L);
    neighbour<2, double>::luaopen(L);
#else
    neighbour<3, float>::luaopen(L);
    neighbour<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class neighbour<3, double>;
template class neighbour<2, double>;
#else
template class neighbour<3, float>;
template class neighbour<2, float>;
#endif

} // namespace host
} // namespace mdsim
} // namespace halmd
