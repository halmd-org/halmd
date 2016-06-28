/*
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

#include <algorithm>
#include <exception>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/max_displacement.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
max_displacement<dimension, float_type>::max_displacement(
    std::shared_ptr<particle_type const> particle
  , std::shared_ptr<box_type const> box
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  // allocate parameters
  , r0_(particle_->nparticle())
  , displacement_(0)
{
}

/**
 * zero maximum displacement
 */
template <int dimension, typename float_type>
void max_displacement<dimension, float_type>::zero()
{
    cache<position_array_type> const& position_cache = particle_->position();
    position_array_type const& position = read_cache(position_cache);
    scoped_timer_type timer(runtime_.zero);
    std::copy(position.begin(), position.end(), r0_.begin());
    displacement_ = 0;
    position_cache_ = position_cache;
}

/**
 * compute maximum displacement
 */
template <int dimension, typename float_type>
float_type max_displacement<dimension, float_type>::compute()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (position_cache != position_cache_) {
        position_array_type const& position = read_cache(position_cache);
        size_type const nparticle = particle_->nparticle();

        scoped_timer_type timer(runtime_.compute);

        float_type rr_max = 0;
        for (typename particle_type::size_type i = 0; i < nparticle; ++i) {
            vector_type r = position[i] - r0_[i];
            box_->reduce_periodic(r);
            rr_max = std::max(rr_max, inner_prod(r, r));
        }
        displacement_ = std::sqrt(rr_max);
        position_cache_ = position_cache;
    }
    return displacement_;
}

template <int dimension, typename float_type>
void max_displacement<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<max_displacement>()
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("zero", &runtime::zero)
                        .def_readonly("compute", &runtime::compute)
                ]
                .def_readonly("runtime", &max_displacement::runtime_)
          , def("max_displacement", &std::make_shared<max_displacement
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<box_type const>
               >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_max_displacement(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    max_displacement<3, double>::luaopen(L);
    max_displacement<2, double>::luaopen(L);
#else
    max_displacement<3, float>::luaopen(L);
    max_displacement<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class max_displacement<3, double>;
template class max_displacement<2, double>;
#else
template class max_displacement<3, float>;
template class max_displacement<2, float>;
#endif

} // namespace host
} // namespace mdsim
} // namespace halmd
