/*
 * Copyright © 2014-2015 Nicolas Höft
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

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_groups/from_region.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <stdexcept>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace particle_groups {

template <typename particle_type>
from_region<particle_type>::from_region(
    std::shared_ptr<particle_type const> particle
  , std::shared_ptr<region_type> region
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , region_(region)
  , logger_(logger)
{
}

template <typename particle_type>
cache<typename from_region<particle_type>::array_type> const&
from_region<particle_type>::ordered()
{
    auto const& selection = region_->selection();
    if (selection != ordered_cache_) {
        auto ordered = make_cache_mutable(ordered_);
        LOG_TRACE("ordered sequence of particle indices");
        ordered->resize(selection->size());
        cuda::copy(selection->begin(), selection->end(), ordered->begin());

        ordered_cache_ = selection;
    }
    return ordered_;
}

template <typename particle_type>
cache<typename from_region<particle_type>::array_type> const&
from_region<particle_type>::unordered()
{
    auto const& selection = region_->selection();
    if (selection != unordered_cache_) {
        auto unordered = make_cache_mutable(unordered_);
        LOG_TRACE("unordered sequence of particle indices");

        unordered->resize(selection->size());
        cuda::copy(selection->begin(), selection->end(), unordered->begin());

        // TODO: is radix sort required here?
        radix_sort(
            unordered->begin()
          , unordered->end()
        );

        unordered_cache_ = selection;
    }
    return unordered_;
}

template <typename particle_type>
cache<typename from_region<particle_type>::size_type> const&
from_region<particle_type>::size()
{
    auto size = make_cache_mutable(size_);
    *size = region_->size();
    return size_;
}

template <typename particle_group_type, typename particle_type>
static void
wrap_to_particle(std::shared_ptr<particle_group_type> self, std::shared_ptr<particle_type> particle_src, std::shared_ptr<particle_type> particle_dst)
{
    particle_group_to_particle(*particle_src, *self, *particle_dst);
}

template <typename particle_type>
void from_region<particle_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("particle_groups")
            [
                class_<from_region, particle_group>()
                    .def("to_particle", &wrap_to_particle<from_region<particle_type>, particle_type>)

              , def("from_region", &std::make_shared<from_region<particle_type>
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<region_type>
                  , std::shared_ptr<logger>
                  >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_particle_groups_from_region(lua_State* L)
{
    from_region<particle<3, float>>::luaopen(L);
    from_region<particle<2, float>>::luaopen(L);
    return 0;
}

// explicit instantiation
template class from_region<particle<3, float>>;
template class from_region<particle<2, float>>;

} // namespace particle_groups
} // namespace gpu
} // namespace mdsim
} // namespace halmd
