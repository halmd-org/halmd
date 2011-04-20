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
#include <lua.hpp>
#include <vector>

#include <halmd/observables/blocking_scheme.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd
{
namespace observables { namespace samples
{

class blocking_scheme_base
{
public:
    typedef observables::blocking_scheme::block_index block_index;

    virtual void push_back(block_index i) = 0;
    virtual void pop_front(block_index i) = 0;
    virtual void clear(block_index i) = 0;
    virtual bool full(block_index i) = 0;
    virtual bool empty(block_index i) = 0;
    virtual std::size_t size(block_index i) = 0;
};

template <typename sample_type>
class blocking_scheme
  : public blocking_scheme_base
{
public:
    typedef blocking_scheme_base _Base;
    typedef typename _Base::block_index block_index;
    typedef boost::circular_buffer<sample_type const> block_type;
    typedef typename block_type::iterator block_iterator;
    typedef typename block_type::const_iterator block_const_iterator;

    static void luaopen(lua_State* L, char const* class_name);

    blocking_scheme(boost::shared_ptr<sample_type const> sample);
    virtual void push_back(block_index i);
    virtual void pop_front(block_index i);
    virtual void clear(block_index i);
    virtual bool full(block_index i);
    virtual bool empty(block_index i);
    virtual std::size_t size(block_index i);

    /**
     * This function is inlined by the correlation function.
     */
    block_type const& block(block_index i) const
    {
        assert(i < blocks_.size());
        return blocks_[i];
    }

private:
    boost::shared_ptr<sample_type const> sample_;
    std::vector<block_type> blocks_;
};

template <typename sample_type>
blocking_scheme<sample_type>::blocking_scheme(boost::shared_ptr<sample_type const> sample)
  : sample_(sample)
{
}

template <typename sample_type>
void blocking_scheme<sample_type>::push_back(block_index i)
{
    assert(i < blocks_.size());
    blocks_[i].push_back(*sample_);
}

template <typename sample_type>
void blocking_scheme<sample_type>::pop_front(block_index i)
{
    assert(i < blocks_.size());
    blocks_[i].pop_front();
}

template <typename sample_type>
void blocking_scheme<sample_type>::clear(block_index i)
{
    assert(i < blocks_.size());
    blocks_[i].clear();
}

template <typename sample_type>
bool blocking_scheme<sample_type>::full(block_index i)
{
    assert(i < blocks_.size());
    return blocks_[i].full();
}

template <typename sample_type>
bool blocking_scheme<sample_type>::empty(block_index i)
{
    assert(i < blocks_.size());
    return blocks_[i].empty();
}

template <typename sample_type>
std::size_t blocking_scheme<sample_type>::size(block_index i)
{
    assert(i < blocks_.size());
    return blocks_[i].size();
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
                    .def(constructor<boost::shared_ptr<sample_type const> >())
            ]
        ]
    ];
}

}} // namespace observables::samples

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SAMPLES_BLOCKING_SCHEME_HPP */
