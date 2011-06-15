/*
 * Copyright © 2011  Peter Colberg and Felix Höfling
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
#include <stdint.h> // uint64_t
#include <vector>

#include <halmd/utility/lua/lua.hpp>

namespace halmd
{
namespace observables { namespace samples
{

/**
 * Abstract base class that defines the interface of the blocking scheme
 * that represents the input data on different coarse-graining levels.
 * Each level roughly models a circular buffer.
 */
class blocking_scheme_base
{
public:
    /** append the current input data to level 'index' */
    virtual void push_back(std::size_t index) = 0;
    /** drop the first entry at level 'index' */
    virtual void pop_front(std::size_t index) = 0;
    /** clear the data of level 'index' */
    virtual void clear(std::size_t index) = 0;
    /** returns true if level 'index' is full */
    virtual bool full(std::size_t index) const = 0;
    /** returns true if level 'index' contains no data */
    virtual bool empty(std::size_t index) const = 0;
    /** returns number of data points stored at level 'index' */
    virtual std::size_t size(std::size_t index) const = 0;
    /** returns number of coarse-graining levels */
    virtual std::size_t count() const = 0;
    /** returns size of coarse-graining blocks */
    virtual std::size_t block_size() const = 0;
    /** returns time stamp (aka integration step) of the current sample */
    virtual uint64_t timestamp() const = 0;
};

/**
 * Represents a set of coarse-grained blocks of input samples
 * of type sample_type, e.g., phase space or density modes.
 */
template <typename sample_type>
class blocking_scheme
  : public blocking_scheme_base
{
public:
    typedef boost::circular_buffer<sample_type> block_type;
    typedef typename block_type::iterator block_iterator;
    typedef typename block_type::const_iterator block_const_iterator;

    static void luaopen(lua_State* L, char const* scope);

    blocking_scheme(
        boost::shared_ptr<sample_type const> sample
      , std::size_t count
      , std::size_t size
    );
    virtual void push_back(std::size_t index);
    virtual void pop_front(std::size_t index);
    virtual void clear(std::size_t index);
    virtual bool full(std::size_t index) const;
    virtual bool empty(std::size_t index) const;
    virtual std::size_t size(std::size_t index) const;
    virtual std::size_t count() const;
    virtual std::size_t block_size() const;
    virtual uint64_t timestamp() const;

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


/**
 * @param sample shared pointer to the current input sample
 * @param count  number of coarse-graining levels
 * @param size   maximum size of each coarse-graining level
 */
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
bool blocking_scheme<sample_type>::full(std::size_t index) const
{
    assert(index < blocks_.size());
    return blocks_[index].full();
}

template <typename sample_type>
bool blocking_scheme<sample_type>::empty(std::size_t index) const
{
    assert(index < blocks_.size());
    return blocks_[index].empty();
}

template <typename sample_type>
std::size_t blocking_scheme<sample_type>::size(std::size_t index) const
{
    assert(index < blocks_.size());
    return blocks_[index].size();
}

template <typename sample_type>
std::size_t blocking_scheme<sample_type>::count() const
{
    return blocks_.size();
}

template <typename sample_type>
std::size_t blocking_scheme<sample_type>::block_size() const
{
    assert(!blocks_.empty());
    return blocks_[0].capacity(); // choose level 0 as representative
}

template <typename sample_type>
uint64_t blocking_scheme<sample_type>::timestamp() const
{
    // return integration step at which sample data were taken
    return sample_->step;
}

template <typename sample_type>
static char const* sample_name_wrapper(blocking_scheme<sample_type> const&)
{
    return sample_type::class_name();
}

template <typename sample_type>
void blocking_scheme<sample_type>::luaopen(lua_State* L, char const* scope)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_(scope)
            [
                namespace_("samples")
                [
                    namespace_("blocking_scheme")
                    [
                        class_<blocking_scheme, boost::shared_ptr<_Base>, _Base>(sample_type::class_name())
                            .def(constructor<
                                boost::shared_ptr<sample_type const>
                              , std::size_t
                              , std::size_t
                            >())
                            .property("sample_name", &sample_name_wrapper<sample_type>)
                    ]
                ]
            ]
        ]
    ];
}

}} // namespace observables::samples

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SAMPLES_BLOCKING_SCHEME_HPP */
