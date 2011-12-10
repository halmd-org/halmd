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

#include <halmd/mdsim/host/neighbours/from_binning.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace neighbours {

/**
 * construct neighbour list module
 *
 * @param particle mdsim::host::particle instance
 * @param box mdsim::box instance
 * @param cutoff force cutoff radius
 * @param skin neighbour list skin
 */
template <int dimension, typename float_type>
from_binning<dimension, float_type>::from_binning(
    shared_ptr<particle_type const> particle1
  , shared_ptr<particle_type const> particle2
  , shared_ptr<binning_type const> binning1
  , shared_ptr<binning_type const> binning2
  , shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , double skin
  , shared_ptr<logger> logger
)
  // dependency injection
  : particle1_(particle1)
  , particle2_(particle2)
  , binning1_(binning1)
  , binning2_(binning2)
  , box_(box)
  , logger_(logger)
  // allocate parameters
  , neighbour_(particle1_->nbox)
  , r_skin_(skin)
  , rr_cut_skin_(particle1_->ntype, particle1_->ntype)
{
    matrix_type r_cut_skin(r_cut.size1(), r_cut.size2());
    typename matrix_type::value_type r_cut_max = 0;
    for (size_t i = 0; i < r_cut.size1(); ++i) {
        for (size_t j = 0; j < r_cut.size2(); ++j) {
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
void from_binning<dimension, float_type>::update()
{
    // Emit on_prepend_update signal, which may be connected e.g. to the
    // binning update slot. We don't call binning::update directly, since
    // the order of calls is setup at the Lua level, and it allows us to
    // pass binning as a const dependency.
    on_prepend_update_();

    LOG_TRACE("update neighbour lists");

    scoped_timer_type timer(runtime_.update);

    cell_size_type const& ncell = binning1_->ncell();
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

    on_append_update_();
}

/**
 * Update neighbour lists for a single cell
 */
template <int dimension, typename float_type>
void from_binning<dimension, float_type>::update_cell_neighbours(cell_size_type const& i)
{
    cell_size_type const& ncell = binning1_->ncell();
    cell_lists const& cell1 = binning1_->cell();
    cell_lists const& cell2 = binning2_->cell();

    BOOST_FOREACH(size_t p, cell1(i)) {
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
                        compute_cell_neighbours<false>(p, cell2(k));
                    }
                }
                else {
                    // visit half of 8 neighbour cells due to pair potential
                    if (j[0] == 0 && j[1] == 0) {
                        goto self;
                    }
                    // update neighbour list of particle
                    cell_size_type k = element_mod(static_cast<cell_size_type>(static_cast<cell_diff_type>(i + ncell) + j), ncell);
                    compute_cell_neighbours<false>(p, cell2(k));
                }
            }
        }
self:
        // visit this cell
        compute_cell_neighbours<true>(p, cell2(i));
    }
}

/**
 * Update neighbour list of particle
 */
template <int dimension, typename float_type>
template <bool same_cell>
void from_binning<dimension, float_type>::compute_cell_neighbours(size_t i, cell_list const& c)
{
    vector<vector_type> const& r1 = particle1_->r;
    vector<vector_type> const& r2 = particle2_->r;
    vector<unsigned int> const& type1 = particle1_->type;
    vector<unsigned int> const& type2 = particle2_->type;

    BOOST_FOREACH(size_t j, c) {
        // skip identical particle and particle pair permutations if same cell
        if (same_cell && particle1_ == particle2_ && j <= i) {
            continue;
        }

        // particle distance vector
        vector_type r = r1[i] - r2[j];
        box_->reduce_periodic(r);
        // particle types
        size_t a = type1[i];
        size_t b = type2[j];
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
void from_binning<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("from_binning_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("neighbours")
                [
                    class_<from_binning, shared_ptr<mdsim::neighbour>, mdsim::neighbour>(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type const>
                          , shared_ptr<particle_type const>
                          , shared_ptr<binning_type const>
                          , shared_ptr<binning_type const>
                          , shared_ptr<box_type const>
                          , matrix_type const&
                          , double
                          , shared_ptr<logger_type>
                         >())
                        .property("r_skin", &from_binning::r_skin)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("update", &runtime::update)
                        ]
                        .def_readonly("runtime", &from_binning::runtime_)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_neighbours_from_binning(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    from_binning<3, double>::luaopen(L);
    from_binning<2, double>::luaopen(L);
#else
    from_binning<3, float>::luaopen(L);
    from_binning<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class from_binning<3, double>;
template class from_binning<2, double>;
#else
template class from_binning<3, float>;
template class from_binning<2, float>;
#endif

} // namespace neighbours
} // namespace host
} // namespace mdsim
} // namespace halmd
