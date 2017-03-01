/*
 * Copyright © 2013-2015  Nicolas Höft
 * Copyright © 2015       Felix Höfling
 * Copyright © 2008-2011  Peter Colberg
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

#include <halmd/mdsim/host/neighbours/from_binning.hpp>
#include <halmd/utility/lua/lua.hpp>

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
    std::pair<std::shared_ptr<particle_type const>, std::shared_ptr<particle_type const>> particle
  , std::pair<std::shared_ptr<binning_type>, std::shared_ptr<binning_type>> binning
  , std::pair<std::shared_ptr<displacement_type>, std::shared_ptr<displacement_type>> displacement
  , std::shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , double skin
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle1_(particle.first)
  , particle2_(particle.second)
  , binning2_(binning.second)
  , displacement1_(displacement.first)
  , displacement2_(displacement.second)
  , box_(box)
  , logger_(logger)
  // allocate parameters
  , neighbour_(particle1_->nparticle())
  , r_skin_(skin)
  , rr_cut_skin_(particle1_->nspecies(), particle2_->nspecies())
{
    matrix_type r_cut_skin(r_cut.size1(), r_cut.size2());
    typename matrix_type::value_type r_cut_max = 0;
    for (size_t i = 0; i < r_cut.size1(); ++i) {
        for (size_t j = 0; j < r_cut.size2(); ++j) {
            r_cut_skin(i, j) = r_cut(i, j) + r_skin_;
            rr_cut_skin_(i, j) = std::pow(r_cut_skin(i, j), 2);
            r_cut_max = std::max(r_cut_skin(i, j), r_cut_max);
        }
    }

    LOG("neighbour list skin: " << r_skin_);
}

template <int dimension, typename float_type>
cache<std::vector<typename from_binning<dimension, float_type>::neighbour_list>> const&
from_binning<dimension, float_type>::lists()
{
    cache<reverse_tag_array_type> const& reverse_tag_cache1 = particle1_->reverse_tag();
    cache<reverse_tag_array_type> const& reverse_tag_cache2 = particle2_->reverse_tag();

    auto current_cache = std::tie(reverse_tag_cache1, reverse_tag_cache2);

    if (neighbour_cache_ != current_cache || displacement1_->compute() > r_skin_ / 2
        || displacement2_->compute() > r_skin_ / 2) {
        on_prepend_update_();
        update();
        displacement1_->zero();
        displacement2_->zero();
        neighbour_cache_ = current_cache;
        on_append_update_();
    }
    return neighbour_;
}

/**
 * Test compatibility of binning parameters with this neighbour list algorithm
 *
 * The binning module is required to have at least 3 cells in each spatial
 * direction in order to be used with the neighbour module.
 */
template <int dimension, typename float_type>
bool from_binning<dimension, float_type>::is_binning_compatible(
    std::shared_ptr<binning_type const> binning1
  , std::shared_ptr<binning_type const> binning2
)
{
    auto ncell = binning2->ncell();
    return *std::min_element(ncell.begin(), ncell.end()) >= 3;
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

    auto const& position1 = read_cache(particle1_->position());
    auto const& cell2 = read_cache(binning2_->cell());
    size_type nparticle1 = particle1_->nparticle();
    cell_size_type const& ncell = binning2_->ncell();

    auto neighbour = make_cache_mutable(neighbour_);

    LOG_TRACE("update neighbour lists");

    scoped_timer_type timer(runtime_.update);

    for (size_type i = 0; i < nparticle1; ++i) {
        // load first particle
        vector_type const& r1 = position1[i];
        cell_size_type idx = binning2_->index(r1);

        // clear particle's neighbour list
        (*neighbour)[i].clear();

        // update neighbour list by iterating over all neighbouring cells
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
                        cell_size_type k = element_mod(static_cast<cell_size_type>(static_cast<cell_diff_type>(idx + ncell) + j), ncell);
                        add_neighbours<false>(i, cell2(k));
                    }
                }
                else {
                    // visit half of 8 neighbour cells due to pair potential
                    if (j[0] == 0 && j[1] == 0) {
                        goto self;
                    }
                    // update neighbour list of particle
                    cell_size_type k = element_mod(static_cast<cell_size_type>(static_cast<cell_diff_type>(idx + ncell) + j), ncell);
                    add_neighbours<false>(i, cell2(k));
                }
            }
        }
self:
        // visit this cell
        add_neighbours<true>(i, cell2(idx));
    }
}

/**
 * Add neighbours of particle i from particles in cell c
 */
template <int dimension, typename float_type>
template <bool same_cell>
void from_binning<dimension, float_type>::add_neighbours(size_t i, cell_list const& c)
{
    auto neighbour = make_cache_mutable(neighbour_);

    position_array_type const& position1 = read_cache(particle1_->position());
    position_array_type const& position2 = read_cache(particle2_->position());
    species_array_type const& species1 = read_cache(particle1_->species());
    species_array_type const& species2 = read_cache(particle2_->species());

    for (size_type j : c) {
        // skip identical particle and particle pair permutations if same cell
        if (same_cell && particle1_ == particle2_ && j <= i) {
            continue;
        }

        // particle distance vector
        vector_type r = position1[i] - position2[j];
        box_->reduce_periodic(r);
        // particle types
        species_type a = species1[i];
        species_type b = species2[j];
        // squared particle distance
        float_type rr = inner_prod(r, r);

        // enforce cutoff radius with neighbour list skin
        if (rr >= rr_cut_skin_(a, b)) {
            continue;
        }

        // add particle to neighbour list
        (*neighbour)[i].push_back(j);
    }
}

template <int dimension, typename float_type>
void from_binning<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("neighbours")
            [
                class_<from_binning, _Base>()
                    .property("r_skin", &from_binning::r_skin)
                    .def("on_prepend_update", &from_binning::on_prepend_update)
                    .def("on_append_update", &from_binning::on_append_update)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("update", &runtime::update)
                    ]
                    .def_readonly("runtime", &from_binning::runtime_)
              , def("from_binning", &std::make_shared<from_binning
                  , std::pair<std::shared_ptr<particle_type const>, std::shared_ptr<particle_type const>>
                  , std::pair<std::shared_ptr<binning_type>, std::shared_ptr<binning_type>>
                  , std::pair<std::shared_ptr<displacement_type>, std::shared_ptr<displacement_type>>
                  , std::shared_ptr<box_type const>
                  , matrix_type const&
                  , double
                  , std::shared_ptr<logger>
                  >)
              , def("is_binning_compatible", &from_binning::is_binning_compatible)
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
