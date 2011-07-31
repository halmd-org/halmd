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

#include <algorithm>
#include <boost/bind.hpp>
#include <exception>

#include <halmd/mdsim/host/binning.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

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
binning<dimension, float_type>::binning(
    shared_ptr<particle_type const> particle
  , shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , float_type skin
  , shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , logger_(logger)
  // allocate parameters
  , r_skin_(skin)
{
    matrix_type r_cut_skin(particle->ntype, particle->ntype);
    typename matrix_type::value_type r_cut_max = 0;
    for (size_t i = 0; i < particle->ntype; ++i) {
        for (size_t j = i; j < particle->ntype; ++j) {
            r_cut_skin(i, j) = r_cut(i, j) + r_skin_;
            r_cut_max = max(r_cut_skin(i, j), r_cut_max);
        }
    }
    vector_type L = box->length();
    ncell_ = static_cast<cell_size_type>(L / r_cut_max);
    if (*min_element(ncell_.begin(), ncell_.end()) < 3) {
        LOG_DEBUG("number of cells per dimension (untruncated): " << L / r_cut_max);
        throw logic_error("less than least 3 cells per dimension");
    }
    cell_.resize(ncell_);
    cell_length_ = element_div(L, static_cast<vector_type>(ncell_));

    LOG("neighbour list skin: " << r_skin_);
    LOG("number of cells per dimension: " << ncell_);
    LOG("cell edge lengths: " << cell_length_);
}

/**
 * Update cell lists
 */
template <int dimension, typename float_type>
void binning<dimension, float_type>::update()
{
    LOG_TRACE("update cell lists");

    scoped_timer_type timer(runtime_.update);

    // empty cell lists without memory reallocation
    for_each(cell_.data(), cell_.data() + cell_.num_elements(), bind(&cell_list::clear, _1));
    // add particles to cells
    for (size_t i = 0; i < particle_->nbox; ++i) {
        vector_type const& r = particle_->r[i];
        cell_size_type index = element_mod(static_cast<cell_size_type>(element_div(r, cell_length_) + static_cast<vector_type>(ncell_)), ncell_);
        cell_(index).push_back(i);
    }
}

template <typename binning_type>
typename signal<void ()>::slot_function_type
wrap_update(shared_ptr<binning_type> binning)
{
    return bind(&binning_type::update, binning);
}

template <int dimension, typename float_type>
void binning<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("binning_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                class_<binning, shared_ptr<binning> >(class_name.c_str())
                    .def(constructor<
                        shared_ptr<particle_type const>
                      , shared_ptr<box_type const>
                      , matrix_type const&
                      , float_type
                      , shared_ptr<logger_type>
                     >())
                    .property("update", &wrap_update<binning>)
                    .property("r_skin", &binning::r_skin)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("update", &runtime::update)
                    ]
                    .def_readonly("runtime", &binning::runtime_)
            ]
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
