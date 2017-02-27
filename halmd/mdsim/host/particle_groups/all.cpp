/*
 * Copyright © 2012 Peter Colberg
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

#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <algorithm>

namespace halmd {
namespace mdsim {
namespace host {
namespace particle_groups {

template <typename particle_type>
all<particle_type>::all(std::shared_ptr<particle_type const> particle)
  : particle_(particle)
  , unordered_(particle_->nparticle())
  , ordered_(particle_->nparticle())
  , size_(particle_->nparticle())
{
    auto unordered = make_cache_mutable(unordered_);
    std::iota(unordered->begin(), unordered->end(), 0);
}

template <typename particle_type>
cache<typename all<particle_type>::array_type> const&
all<particle_type>::ordered()
{
    if (!(ordered_observer_ == particle_->reverse_id())) {
        get_reverse_id(*particle_, make_cache_mutable(ordered_)->begin());
        ordered_observer_ = particle_->reverse_id();
    }
    return ordered_;
}

template <typename particle_type>
cache<typename all<particle_type>::array_type> const&
all<particle_type>::unordered()
{
    return unordered_;
}

template <typename particle_type>
cache<typename all<particle_type>::size_type> const&
all<particle_type>::size()
{
    return size_;
}

template <typename particle_type>
void all<particle_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("particle_groups")
            [
                class_<all, particle_group>()

              , def("all", &std::make_shared<all, std::shared_ptr<particle_type const>>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_particle_groups_all(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    all<particle<3, double>>::luaopen(L);
    all<particle<2, double>>::luaopen(L);
#else
    all<particle<3, float>>::luaopen(L);
    all<particle<2, float>>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class all<particle<3, double>>;
template class all<particle<2, double>>;
#else
template class all<particle<3, float>>;
template class all<particle<2, float>>;
#endif

} // namespace particle_groups
} // namespace host
} // namespace mdsim
} // namespace halmd
