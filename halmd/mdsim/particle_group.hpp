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

#include <boost/make_shared.hpp>
#include <cstddef>
#include <lua.hpp>
#include <utility>

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
     * Returns underlying particle instance.
     */
    virtual boost::shared_ptr<particle_type const> particle() const = 0;

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);
};

template <typename particle_type>
void particle_group<particle_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L)
    [
        class_<particle_group>()
            .property("particle", &particle_group::particle)
            .property("size", &particle_group::size)
    ];
}

/**
 * Select particles of a given particle instance by a contiguous range of particle tags.
 */
template <typename particle_type>
class particle_group_from_range
  : public particle_group<particle_type>
{
private:
    typedef particle_group<particle_type> _Base;

public:
    typedef typename _Base::iterator iterator;

    /**
     * Select by tag range [begin, end).
     */
    particle_group_from_range(
        boost::shared_ptr<particle_type const> particle
      , std::size_t begin
      , std::size_t end
    );

    /**
     * Returns iterator to first element of array mapping particle tags to array.
     */
    virtual iterator begin() const;

    /**
     * Returns iterator past last element of array mapping particle tags to array.
     */
    virtual iterator end() const;

    /**
     * Returns number of particles.
     */
    virtual std::size_t size() const;

    /**
     * Returns underlying particle instance.
     */
    virtual boost::shared_ptr<particle_type const> particle() const;

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** particle instance */
    boost::shared_ptr<particle_type const> particle_;
    /** first tag of range */
    std::size_t begin_;
    /** past last tag of range */
    std::size_t end_;
};

template <typename particle_type>
particle_group_from_range<particle_type>::particle_group_from_range(
    boost::shared_ptr<particle_type const> particle
  , std::size_t begin
  , std::size_t end
)
  : particle_(particle)
  , begin_(begin)
  , end_(end)
{
    if (end_ < begin_) {
        throw std::logic_error("particle_group: inverse tag ranges not allowed");
    }
    if (end_ > particle_->reverse_tag().size()) {
        throw std::logic_error("particle_group: tag range exceeds particle array");
    }
}

template <typename particle_type>
typename particle_group_from_range<particle_type>::iterator
particle_group_from_range<particle_type>::begin() const
{
    return particle_->reverse_tag().begin() + begin_;
}

template <typename particle_type>
typename particle_group_from_range<particle_type>::iterator
particle_group_from_range<particle_type>::end() const
{
    return particle_->reverse_tag().begin() + end_;
}

template <typename particle_type>
std::size_t
particle_group_from_range<particle_type>::size() const
{
    return end_ - begin_;
}

template <typename particle_type>
boost::shared_ptr<particle_type const>
particle_group_from_range<particle_type>::particle() const
{
    return particle_;
}

/**
 * Convert from range of 1-based tags to 0-based iterator range
 */
template <typename particle_type>
static boost::shared_ptr<particle_group_from_range<particle_type> >
wrap_particle_group_from_range(
    boost::shared_ptr<particle_type const> particle
  , std::size_t first
  , std::size_t last
)
{
    return boost::make_shared<particle_group_from_range<particle_type> >(
        particle
      , first - 1
      , last
    );
}

template <typename particle_type>
void particle_group_from_range<particle_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<particle_group_from_range, _Base>()

          , def("particle_group_from_range", &wrap_particle_group_from_range<particle_type>)
        ]
    ];
}

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_PARTICLE_GROUP_HPP */
