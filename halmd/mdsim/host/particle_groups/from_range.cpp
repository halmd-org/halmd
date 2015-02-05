/*
 * Copyright © 2012 Peter Colberg
 * Copyright © 2012 Felix Höfling
 * Copyright © 2015 Nicolas Höft
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

#include <halmd/algorithm/host/radix_sort.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/from_range.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <algorithm>
#include <stdexcept>

namespace halmd {
namespace mdsim {
namespace host {
namespace particle_groups {

template <typename particle_type>
from_range<particle_type>::from_range(
    std::shared_ptr<particle_type const> particle
  , range_type const& range
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , range_(check_range(range))
  , logger_(logger)
  , ordered_(range.second - range.first)
  , unordered_(range.second - range.first)
  , size_(range.second - range.first)
{
}

template <typename particle_type>
typename from_range<particle_type>::range_type const&
from_range<particle_type>::check_range(range_type const& range)
{
    if (range.second <= range.first) {
        throw std::invalid_argument("particle_group: inverse tag ranges not allowed");
    }
    if (range.second > particle_->nparticle()) {
        throw std::invalid_argument("particle_group: tag range exceeds particle array");
    }
    return range;
}

template <typename particle_type>
cache<typename from_range<particle_type>::array_type> const&
from_range<particle_type>::ordered()
{
    cache<array_type> const& reverse_tag_cache = particle_->reverse_tag();
    if (ordered_cache_ != reverse_tag_cache) {
        LOG_TRACE("ordered sequence of particle indices");
        array_type const& reverse_tag = read_cache(reverse_tag_cache);
        auto ordered = make_cache_mutable(ordered_);
        std::copy(
            reverse_tag.begin() + range_.first
          , reverse_tag.begin() + range_.second
          , ordered->begin()
        );
        ordered_cache_ = reverse_tag_cache;
    }
    return ordered_;
}

template <typename particle_type>
cache<typename from_range<particle_type>::array_type> const&
from_range<particle_type>::unordered()
{
    cache<array_type> const& reverse_tag_cache = particle_->reverse_tag();
    if (unordered_cache_ != reverse_tag_cache) {
        LOG_TRACE("unordered sequence of particle indices");
        array_type const& reverse_tag = read_cache(reverse_tag_cache);
        auto unordered = make_cache_mutable(unordered_);
        std::copy(
            reverse_tag.begin() + range_.first
          , reverse_tag.begin() + range_.second
          , unordered->begin()
        );
        radix_sort(
            unordered->begin()
          , unordered->end()
        );
        unordered_cache_ = reverse_tag_cache;
    }
    return unordered_;
}

template <typename particle_type>
cache<typename from_range<particle_type>::size_type> const&
from_range<particle_type>::size()
{
    return size_;
}

/**
 * This function serves as a Lua wrapper around the C++ constructor,
 * converting from a 1-based particle tag range to a 0-based range.
 */
template <typename particle_type>
static std::shared_ptr<from_range<particle_type> >
wrap_from_range(
    std::shared_ptr<particle_type const> particle
  , typename from_range<particle_type>::range_type const& range
  , std::shared_ptr<logger> logger
)
{
    return std::make_shared<from_range<particle_type> >(
        particle
      , std::make_pair(range.first - 1, range.second)
      , logger
    );
}

template <typename particle_group_type, typename particle_type>
static void
wrap_to_particle(std::shared_ptr<particle_group_type> self, std::shared_ptr<particle_type> particle_src, std::shared_ptr<particle_type> particle_dst)
{
    particle_group_to_particle(*particle_src, *self, *particle_dst);
}


template <typename particle_type>
void from_range<particle_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("particle_groups")
            [
                class_<from_range, particle_group>()
                   .def("to_particle", &wrap_to_particle<from_range<particle_type>, particle_type>)

              , def("from_range", &wrap_from_range<particle_type>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_particle_groups_from_range(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    from_range<particle<3, double>>::luaopen(L);
    from_range<particle<2, double>>::luaopen(L);
#else
    from_range<particle<3, float>>::luaopen(L);
    from_range<particle<2, float>>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class from_range<particle<3, double>>;
template class from_range<particle<2, double>>;
#else
template class from_range<particle<3, float>>;
template class from_range<particle<2, float>>;
#endif

} // namespace particle_groups
} // namespace host
} // namespace mdsim
} // namespace halmd
