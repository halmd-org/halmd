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

#include <halmd/algorithm/gpu/iota.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_groups/all.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace particle_groups {

template <typename particle_type>
all<particle_type>::all(std::shared_ptr<particle_type const> particle)
  : particle_(particle)
  , unordered_(particle_->nparticle())
  , ordered_(particle_->nparticle())
  , size_(particle_->nparticle())
{
    auto unordered = make_cache_mutable(unordered_);
    halmd::iota(unordered->begin(), unordered->end(), 0);
}

template <typename particle_type>
cache<typename all<particle_type>::array_type> const&
all<particle_type>::ordered()
{
    if (!(ordered_observer_ == particle_->reverse_id())) {
        auto reverse_id = read_cache(particle_->reverse_id());
        cuda::copy(reverse_id.begin(), reverse_id.begin() + particle_->nparticle(), make_cache_mutable(ordered_)->begin());
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

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_particle_groups_all(lua_State* L)
{
    all<particle<3, float>>::luaopen(L);
    all<particle<2, float>>::luaopen(L);
    return 0;
}

// explicit instantiation
template class all<particle<3, float>>;
template class all<particle<2, float>>;

} // namespace particle_groups
} // namespace gpu
} // namespace mdsim
} // namespace halmd
