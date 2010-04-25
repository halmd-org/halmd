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
#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/options.hpp>

namespace halmd
{
namespace utility { namespace detail
{

template <typename T, typename Enable = void>
class builder
  : public builder<typename T::_Base>
{
public:
    typedef builder<typename T::_Base> _Base;
    typedef typename _Base::resolve_ptr resolve_ptr;
    typedef typename _Base::create_ptr create_ptr;
    typedef typename _Base::base_ptr base_ptr;
    typedef typename _Base::builder_base_ptr builder_base_ptr;

    builder()
      : _Base(&T::resolve, &create_, priority_)
    {}

    template <typename builder_type>
    builder(builder_type const& builder_, typename boost::enable_if<
      boost::is_base_of<builder_type, builder>, void>::type* dummy = 0)
      : _Base(builder_)
    {}

    virtual ~builder() {}

    virtual bool is_base_of(builder_base_ptr const& builder_)
    {
        return boost::dynamic_pointer_cast<builder>(builder_);
    }

protected:
    builder(resolve_ptr resolve, create_ptr create, int priority)
      : _Base(resolve, create, priority)
    {}

    enum { priority_ = _Base::priority_ + 1 };

private:
    static base_ptr create_(options const& vm)
    {
        return base_ptr(new T(vm));
    }
};

template <typename T>
class builder<T, typename boost::enable_if<
  boost::is_same<T, typename T::module_type> >::type>
{
public:
    typedef void (*resolve_ptr)(options const&);
    typedef typename T::module_type base_type;
    typedef boost::shared_ptr<base_type> base_ptr;
    typedef base_ptr (*create_ptr)(options const&);
    typedef boost::shared_ptr<builder> builder_base_ptr;

    resolve_ptr const resolve;
    create_ptr const create;
    int const priority;

    builder()
      : resolve(&T::resolve)
      , create(&create_)
      , priority(priority_)
    {}

    template <typename builder_type>
    builder(builder_type const& builder_, typename boost::enable_if<
      boost::is_base_of<builder_type, builder>, void>::type* dummy = 0)
      : resolve(builder_.resolve)
      , create(builder_.create)
      , priority(builder_.priority)
    {}

    virtual ~builder() {}

    virtual bool is_base_of(builder_base_ptr const& builder_)
    {
        return boost::dynamic_pointer_cast<builder>(builder_);
    }

    struct less
    {
        bool operator()(builder_base_ptr const& a, builder_base_ptr const& b) const
        {
            if (a->priority == b->priority) {
                return a->resolve != b->resolve;
            }
            return a->priority > b->priority;
        }
    };

protected:
    builder(resolve_ptr resolve, create_ptr create, int priority)
      : resolve(resolve)
      , create(create)
      , priority(priority)
    {}

    enum { priority_ = 0 };

private:
    static base_ptr create_(options const& vm)
    {
        return base_ptr(new T(vm));
    }
};

}} // namespace utility::detail

} // namespace halmd

#endif /* ! HALMD_UTILITY_DETAIL_BUILDER_HPP */
