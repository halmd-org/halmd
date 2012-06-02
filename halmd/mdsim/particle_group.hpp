/*
 * Copyright © 2012 Felix Höfling
 * Copyright © 2012 Peter Colberg
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

#ifndef HALMD_MDSIM_PARTICLE_GROUP_HPP
#define HALMD_MDSIM_PARTICLE_GROUP_HPP

#include <cstddef>
#include <lua.hpp>

#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {

/**
 * A particle group represents a subset of particles, which is defined
 * by an instance of particle together with a sequence of indices.
 */
template <typename particle_type>
class particle_group
{
public:
    typedef typename particle_type::reverse_tag_array_type::const_iterator iterator;

    /**
     * Returns iterator to first element of index sequence.
     */
    virtual iterator begin() const = 0;

    /**
     * Returns iterator one past last element of index sequence.
     */
    virtual iterator end() const = 0;

    /**
     * Returns number of particles.
     */
    virtual std::size_t size() const = 0;

    /**
     * Returns reference to underlying particle instance.
     */
    virtual particle_type& particle() = 0;

    /**
     * Returns const reference to underlying particle instance.
     */
    virtual particle_type const& particle() const = 0;

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);
};

template <typename particle_type>
static particle_type&
wrap_particle(particle_group<particle_type>& self)
{
    return self.particle();
}

template <typename particle_type>
void particle_group<particle_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L)
    [
        class_<particle_group>()
            .property("particle", &wrap_particle<particle_type>)
            .property("size", &particle_group::size)
    ];
}

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_PARTICLE_GROUP_HPP */
