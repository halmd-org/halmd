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

#ifndef HALMD_UTILITY_MODULE_BUILDER_HPP
#define HALMD_UTILITY_MODULE_BUILDER_HPP

#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_object.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/utility/options.hpp>

namespace halmd
{
namespace utility { namespace module
{

// import into namespace
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::enable_if;
using boost::is_object;

/**
 * The builder serves as an abstract base hierarchy to the module
 * wrapper. The builder has three variants:
 */
template <typename T = void, typename Enable = void>
struct builder;

/**
 * A type-independent base class, which is used in the factory to
 * store pointers of arbitrary module wrapper instances in a
 * container.
 */
template <>
struct builder<>
{
    virtual ~builder() {}
    virtual void options(po::options_description& desc) = 0;
    virtual void resolve(po::options const& vm) = 0;
    virtual std::string name() = 0;
    virtual std::type_info const& type() = 0;
};

/**
 * A base template, which corresponds to the base class of a
 * module (e.g. force, particle). The module fetch function
 * casts type-independent builder pointers to this type to
 * create module instances.
 */
template <typename T, typename Enable>
struct builder
  : public builder<>
{
    typedef T _Module_base;
    virtual shared_ptr<_Module_base> fetch(po::options const& vm) = 0;

    /**
     * assemble module options
     */
    virtual void options(po::options_description& desc)
    {
        T::options(desc);
    }

    /**
     * resolve class dependencies
     */
    virtual void resolve(po::options const& vm)
    {
        T::resolve(vm);
    }
};

/**
 * A deriving template to the base template.
 */
template <typename T>
struct builder<T, typename enable_if<is_object<typename T::_Base> >::type>
  : public builder<typename T::_Base>
{
    typedef typename T::_Base _Base;

    /**
     * assemble module options
     */
    virtual void options(po::options_description& desc)
    {
        builder<_Base>::options(desc);
        // check for inherited base module function
        if (T::options != _Base::options) {
            T::options(desc);
        }
    }

    /**
     * resolve class dependencies
     */
    virtual void resolve(po::options const& vm)
    {
        builder<_Base>::resolve(vm);
        // check for inherited base module function
        if (T::resolve != _Base::resolve) {
            T::resolve(vm);
        }
    }
};

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_BUILDER_HPP */
