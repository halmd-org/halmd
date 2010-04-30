/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_UTILITY_DETAIL_BUILDER_HPP
#define HALMD_UTILITY_DETAIL_BUILDER_HPP

#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/utility/options.hpp>
#include <halmd/util/logger.hpp>

namespace halmd
{
namespace utility { namespace detail
{

// import into namespace
using boost::dynamic_pointer_cast;
using boost::shared_ptr;

template <typename T = void, typename Enable = void>
struct builder;

/**
 * Type-independent module base class
 */
template <>
struct builder<>
{
    virtual ~builder() {}
    virtual void options(po::options_description& desc) = 0;
    virtual void resolve(po::options const& vm) = 0;
};

/**
 * Abstract module specification
 */
template <typename T, typename Enable>
struct builder
  : builder<typename T::_Base>
{};

template <typename T>
struct builder<T, typename boost::enable_if<
  boost::is_same<T, typename T::module_type> >::type>
  : public builder<>
{
    typedef typename T::module_type _Base;
    virtual bool rank(shared_ptr<builder<_Base> > const& other) const = 0;
    virtual shared_ptr<_Base> fetch(po::options const& vm) = 0;
};

/**
 * Helper class for builder ordering in STL set
 */
template <typename _Base>
struct _builder_order
{
    typedef shared_ptr<builder<_Base> > pointer;
    bool operator()(pointer const& first, pointer const& second) const
    {
        return first->rank(second);
    }
};

}} // namespace utility::detail

} // namespace halmd

#endif /* ! HALMD_UTILITY_DETAIL_BUILDER_HPP */
