/*
 * Copyright © 2011  Peter Colberg
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

#ifndef HALMD_OBSERVABLES_SAMPLES_BLOCKING_SCHEME_HPP
#define HALMD_OBSERVABLES_SAMPLES_BLOCKING_SCHEME_HPP

// Boost 1.37.0, or patch from http://svn.boost.org/trac/boost/ticket/1852
#include <boost/circular_buffer.hpp>
#include <cassert>
#include <cstddef> // std::size_t
#include <lua.hpp>
#include <vector>

#include <halmd/utility/lua/lua.hpp>

namespace halmd
{
namespace observables { namespace samples
{

class blocking_scheme_base
{
public:
    virtual void push_back(std::size_t index) = 0;
    virtual void pop_front(std::size_t index) = 0;
    virtual void clear(std::size_t index) = 0;
    virtual bool full(std::size_t index) = 0;
    virtual bool empty(std::size_t index) = 0;
    virtual std::size_t size(std::size_t index) = 0;
};

template <typename sample_type>
class blocking_scheme
  : public blocking_scheme_base
{
public:
    typedef boost::circular_buffer<sample_type> block_type;
    typedef typename block_type::iterator block_iterator;
    typedef typename block_type::const_iterator block_const_iterator;

    static void luaopen(lua_State* L, char const* class_name);

    blocking_scheme(
        boost::shared_ptr<sample_type const> sample
      , std::size_t count
      , std::size_t size
    );
    virtual void push_back(std::size_t index);
    virtual void pop_front(std::size_t index);
    virtual void clear(std::size_t index);
    virtual bool full(std::size_t index);
    virtual bool empty(std::size_t index);
    virtual std::size_t size(std::size_t index);

    /**
     * This function is inlined by the correlation function.
     */
    block_type const& index(std::size_t index) const
    {
        assert(index < blocks_.size());
        return blocks_[index];
    }

private:
    typedef blocking_scheme_base _Base;

    boost::shared_ptr<sample_type const> sample_;
    std::vector<block_type> blocks_;
};

template <typename sample_type>
blocking_scheme<sample_type>::blocking_scheme(
    boost::shared_ptr<sample_type const> sample
  , std::size_t count
  , std::size_t size
)
  : sample_(sample)
  , blocks_(count, block_type(size))
{
}

template <typename sample_type>
void blocking_scheme<sample_type>::push_back(std::size_t index)
{
    assert(index < blocks_.size());
    blocks_[index].push_back(*sample_);
}

template <typename sample_type>
void blocking_scheme<sample_type>::pop_front(std::size_t index)
{
    assert(index < blocks_.size());
    blocks_[index].pop_front();
}

template <typename sample_type>
void blocking_scheme<sample_type>::clear(std::size_t index)
{
    assert(index < blocks_.size());
    blocks_[index].clear();
}

template <typename sample_type>
bool blocking_scheme<sample_type>::full(std::size_t index)
{
    assert(index < blocks_.size());
    return blocks_[index].full();
}

template <typename sample_type>
bool blocking_scheme<sample_type>::empty(std::size_t index)
{
    assert(index < blocks_.size());
    return blocks_[index].empty();
}

template <typename sample_type>
std::size_t blocking_scheme<sample_type>::size(std::size_t index)
{
    assert(index < blocks_.size());
    return blocks_[index].size();
}

template <typename sample_type>
void blocking_scheme<sample_type>::luaopen(lua_State* L, char const* class_name)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("samples")
            [
                class_<blocking_scheme, boost::shared_ptr<_Base>, _Base>(class_name)
                    .def(constructor<
                        boost::shared_ptr<sample_type const>
                      , std::size_t
                      , std::size_t
                    >())
            ]
        ]
    ];
}

}} // namespace observables::samples

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SAMPLES_BLOCKING_SCHEME_HPP */
