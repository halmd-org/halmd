/*
 * Copyright © 2008-2011  Peter Colberg
 * Copyright © 2013       Nicolas Höft
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
#include <exception>

#include <halmd/mdsim/host/binning.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

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
binning<dimension, float_type>::binning(
    std::shared_ptr<particle_type const> particle
  , std::shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , float_type skin
  , std::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , logger_(logger)
  // allocate parameters
  , r_skin_(skin)
{
    matrix_type r_cut_skin(r_cut.size1(), r_cut.size2());
    typename matrix_type::value_type r_cut_max = 0;
    for (size_t i = 0; i < r_cut.size1(); ++i) {
        for (size_t j = 0; j < r_cut.size2(); ++j) {
            r_cut_skin(i, j) = r_cut(i, j) + r_skin_;
            r_cut_max = std::max(r_cut_skin(i, j), r_cut_max);
        }
    }
    vector_type L = box->length();
    ncell_ = element_max(static_cast<cell_size_type>(L / r_cut_max), cell_size_type(1));

    auto cell = make_cache_mutable(cell_);
    cell->resize(ncell_);
    cell_length_ = element_div(L, static_cast<vector_type>(ncell_));

    LOG("neighbour list skin: " << r_skin_);
    LOG("number of cells per dimension: " << ncell_);
    LOG("cell edge lengths: " << cell_length_);
}

template <int dimension, typename float_type>
cache<typename binning<dimension, float_type>::array_type> const&
binning<dimension, float_type>::cell()
{
    cache<position_array_type> const& position_cache = particle_->position();
    if (cell_cache_ != position_cache) {
        update();
        cell_cache_ = position_cache;
    }
    return cell_;
}

/**
 * Update cell lists
 */
template <int dimension, typename float_type>
void binning<dimension, float_type>::update()
{
    position_array_type const& position = read_cache(particle_->position());
    size_type const nparticle = particle_->nparticle();

    auto cell = make_cache_mutable(cell_);

    LOG_TRACE("update cell lists");

    scoped_timer_type timer(runtime_.update);

    // empty cell lists without memory reallocation
    std::for_each(
        cell->data()
      , cell->data() + cell->num_elements()
      , [](cell_list& cell) {
            cell.clear();
        }
    );
    // add particles to cells
    for (size_type i = 0; i < nparticle; ++i) {
        vector_type const& r = position[i];
        cell_size_type index = element_mod(static_cast<cell_size_type>(element_div(r, cell_length_) + static_cast<vector_type>(ncell_)), ncell_);
        (*cell)(index).push_back(i);
    }
}

template <int dimension, typename float_type>
void binning<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<binning>()
                .property("r_skin", &binning::r_skin)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("update", &runtime::update)
                ]
                .def_readonly("runtime", &binning::runtime_)
          , def("binning", &std::make_shared<binning
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<box_type const>
                  , matrix_type const&
                  , float_type
                  , std::shared_ptr<logger_type>
              >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_binning(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    binning<3, double>::luaopen(L);
    binning<2, double>::luaopen(L);
#else
    binning<3, float>::luaopen(L);
    binning<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class binning<3, double>;
template class binning<2, double>;
#else
template class binning<3, float>;
template class binning<2, float>;
#endif

} // namespace host
} // namespace mdsim
} // namespace halmd
