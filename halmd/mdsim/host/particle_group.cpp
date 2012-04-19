/*
 * Copyright © 2012  Felix Höfling
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

#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <string>

#include <halmd/mdsim/host/particle_group.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
static int wrap_dimension(particle_group<dimension, float_type> const&)
{
    return dimension;
}

template <int dimension, typename float_type>
void particle_group<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name("particle_group_host_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<particle_group, shared_ptr<particle_group> >(class_name.c_str())
                .property("particle", &particle_group::particle)
                .property("size", &particle_group::size)
                .property("empty", &particle_group::empty)
                .property("dimension", &wrap_dimension<dimension, float_type>)
        ]
    ];
}

template <int dimension, typename float_type>
void particle_group_all<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name("particle_group_all_host_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<particle_group_all, shared_ptr<_Base>, _Base>(class_name.c_str())
          , def("particle_group_all", &make_shared<particle_group_all
              , shared_ptr<particle_type const>
            >)
        ]
    ];
}

template <int dimension, typename float_type>
particle_group_all<dimension, float_type>::particle_group_all(
    shared_ptr<particle_type const> particle
)
  : particle_(particle) {}

template <int dimension, typename float_type>
boost::shared_ptr<typename particle_group_all<dimension, float_type>::particle_type const>
particle_group_all<dimension, float_type>::particle() const
{
    return particle_;
}

template <int dimension, typename float_type>
typename particle_group_all<dimension, float_type>::map_iterator
particle_group_all<dimension, float_type>::map() const
{
    return particle_->reverse_tag().begin();
}

template <int dimension, typename float_type>
unsigned int particle_group_all<dimension, float_type>::size() const
{
    return particle_->nparticle();
}

template <int dimension, typename float_type>
bool particle_group_all<dimension, float_type>::all() const
{
    return true;
}

template <int dimension, typename float_type>
particle_group_from_range<dimension, float_type>::particle_group_from_range(
    shared_ptr<particle_type const> particle
  , unsigned int begin, unsigned int end
)
  : particle_(particle)
  , begin_(begin)
  , end_(end)
{
    if (end_ < begin_) {
        throw std::logic_error("particle_group: inverse tag ranges not allowed.");
    }

    if (end_ > particle_->nparticle()) {
        throw std::logic_error("particle_group: tag range exceeds particle array.");
    }
}

template <int dimension, typename float_type>
shared_ptr<typename particle_group_from_range<dimension, float_type>::particle_type const>
particle_group_from_range<dimension, float_type>::particle() const
{
    return particle_;
}

template <int dimension, typename float_type>
typename particle_group_from_range<dimension, float_type>::map_iterator
particle_group_from_range<dimension, float_type>::map() const
{
    return particle_->reverse_tag().begin() + begin_;
}

template <int dimension, typename float_type>
unsigned int particle_group_from_range<dimension, float_type>::size() const
{
    return end_ - begin_;
}

template <int dimension, typename float_type>
bool particle_group_from_range<dimension, float_type>::all() const
{
    return size() == particle_->nparticle();
}

template <int dimension, typename float_type>
void particle_group_from_range<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name("particle_group_from_range_host_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<particle_group_from_range, shared_ptr<_Base>, _Base>(class_name.c_str())
          , def("particle_group_from_range", &make_shared<particle_group_from_range
              , shared_ptr<particle_type const>
              , unsigned int, unsigned int
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_particle_group(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    particle_group<3, double>::luaopen(L);
    particle_group<2, double>::luaopen(L);
    particle_group_all<3, double>::luaopen(L);
    particle_group_all<2, double>::luaopen(L);
    particle_group_from_range<3, double>::luaopen(L);
    particle_group_from_range<2, double>::luaopen(L);
#else
    particle_group<3, float>::luaopen(L);
    particle_group<2, float>::luaopen(L);
    particle_group_all<3, float>::luaopen(L);
    particle_group_all<2, float>::luaopen(L);
    particle_group_from_range<3, float>::luaopen(L);
    particle_group_from_range<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class particle_group_all<3, double>;
template class particle_group_all<2, double>;
template class particle_group_from_range<3, double>;
template class particle_group_from_range<2, double>;
#else
template class particle_group_all<3, float>;
template class particle_group_all<2, float>;
template class particle_group_from_range<3, float>;
template class particle_group_from_range<2, float>;
#endif

} // namespace host
} // namespace mdsim
} // namespace halmd
