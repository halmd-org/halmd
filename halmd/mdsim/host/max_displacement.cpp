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
{
}

/**
 * zero maximum displacement
 */
template <int dimension, typename float_type>
void max_displacement<dimension, float_type>::zero()
{
    position_array_type const& position = read_cache(particle_->position());
    scoped_timer_type timer(runtime_.zero);
    std::copy(position.begin(), position.end(), r0_.begin());
}

/**
 * compute maximum displacement
 */
template <int dimension, typename float_type>
float_type max_displacement<dimension, float_type>::compute()
{
    position_array_type const& position = read_cache(particle_->position());
    size_type const nparticle = particle_->nparticle();

    scoped_timer_type timer(runtime_.compute);

    float_type rr_max = 0;
    for (typename particle_type::size_type i = 0; i < nparticle; ++i) {
        vector_type r = position[i] - r0_[i];
        box_->reduce_periodic(r);
        rr_max = std::max(rr_max, inner_prod(r, r));
    }
    return std::sqrt(rr_max);
}

template <int dimension, typename float_type>
void max_displacement<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static std::string const class_name("max_displacement_" + std::to_string(dimension));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                class_<max_displacement, std::shared_ptr<max_displacement> >(class_name.c_str())
                    .def(constructor<
                         std::shared_ptr<particle_type const>
                       , std::shared_ptr<box_type const>
                     >())
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("zero", &runtime::zero)
                            .def_readonly("compute", &runtime::compute)
                    ]
                    .def_readonly("runtime", &max_displacement::runtime_)
            ]
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
