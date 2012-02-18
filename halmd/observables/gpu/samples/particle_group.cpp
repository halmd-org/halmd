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
#include <string>

#include <halmd/observables/gpu/samples/particle_group.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace gpu {
namespace samples {

template <int dimension, typename float_type>
particle_group_from_range<dimension, float_type>::particle_group_from_range(
    shared_ptr<particle_type /* FIXME const */> particle
  , unsigned int begin, unsigned int end
)
  : particle_(particle)
  , begin_(begin)
  , end_(end)
{
    if (end_ < begin_) {
        throw std::logic_error("particle_group: inverse tag ranges not allowed.");
    }

    if (end_ > particle_->g_reverse_tag.size()) {
        throw std::logic_error("particle_group: tag range exceeds particle array.");
    }
}

template <int dimension, typename float_type>
unsigned int const* particle_group_all<dimension, float_type>::h_map()
{
    // FIXME using particle_->h_reverse_tag prohibits particle_ to be const
    try {
        cuda::copy(particle_->g_reverse_tag, particle_->h_reverse_tag);
    }
    catch (cuda::error const&) {
        throw;
    }

    return particle_->h_reverse_tag.data();
}

template <int dimension, typename float_type>
unsigned int const* particle_group_from_range<dimension, float_type>::h_map()
{
    // FIXME using particle_->h_reverse_tag prohibits particle_ to be const
    try {
        cuda::copy(particle_->g_reverse_tag, particle_->h_reverse_tag);
    }
    catch (cuda::error const&) {
        throw;
    }

    return particle_->h_reverse_tag.data() + begin_;
}

template <int dimension, typename float_type>
void particle_group<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name("particle_group" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                namespace_("samples")
                [
                    class_<particle_group, shared_ptr<particle_group> >(class_name.c_str())
                        .property("particle", &particle_group::particle)
                        .property("size", &particle_group::size)
                        .property("empty", &particle_group::empty)
                ]
            ]
        ]
    ];
}

template <int dimension, typename float_type>
void particle_group_all<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name("particle_group_all" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                namespace_("samples")
                [
                    class_<particle_group_all, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type /* FIXME const*/>
                        >())
                ]
            ]
        ]
    ];
}

template <int dimension, typename float_type>
void particle_group_from_range<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name("particle_group_from_range" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                namespace_("samples")
                [
                    class_<particle_group_from_range, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type /* FIXME const*/>
                          , unsigned int, unsigned int
                        >())
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_samples_particle_group(lua_State* L)
{
    particle_group<3, float>::luaopen(L);
    particle_group<2, float>::luaopen(L);
    return 0;
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_samples_particle_group_all(lua_State* L)
{
    particle_group_all<3, float>::luaopen(L);
    particle_group_all<2, float>::luaopen(L);
    return 0;
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_samples_particle_group_from_range(lua_State* L)
{
    particle_group_from_range<3, float>::luaopen(L);
    particle_group_from_range<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class particle_group_all<3, float>;
template class particle_group_all<2, float>;
template class particle_group_from_range<3, float>;
template class particle_group_from_range<2, float>;

} // namespace samples
} // namespace gpu
} // namespace observables
} // namespace halmd
