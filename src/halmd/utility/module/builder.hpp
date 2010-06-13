/*
 * Copyright © 2010  Peter Colberg
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
#include <boost/weak_ptr.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/utility/module/demangle.hpp>
#include <halmd/utility/module/exception.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace utility { namespace module
{

// import into namespace
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::weak_ptr;
using boost::enable_if;
using boost::is_object;

/**
 * The builder serves as an abstract base hierarchy to the module
 * wrapper. The builder has three variants:
 */

/**
 * A type-independent base class, which is used in the factory to
 * store pointers of arbitrary module wrapper instances in a
 * container.
 */
struct untyped_builder_base
{
    virtual ~untyped_builder_base() {}
    virtual void options(po::options_description& desc) = 0;
    virtual void resolve(po::options const& vm) = 0;
    virtual std::string name() = 0;
    virtual std::type_info const& type() = 0;
    shared_ptr<po::options> vm;
};

typedef shared_ptr<untyped_builder_base> builder;

/**
 * A base template, which corresponds to the base class of a
 * module (e.g. force, particle). The module fetch function
 * casts type-independent builder pointers to this type to
 * create module instances.
 */
template <typename T, typename Enable = void>
struct typed_builder_base
  : public untyped_builder_base
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
struct typed_builder_base<T, typename enable_if<is_object<typename T::_Base> >::type>
  : public typed_builder_base<typename T::_Base>
{
    typedef typename T::_Base _Base;

    /**
     * assemble module options
     */
    virtual void options(po::options_description& desc)
    {
        typed_builder_base<_Base>::options(desc);
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
        typed_builder_base<_Base>::resolve(vm);
        // check for inherited base module function
        if (T::resolve != _Base::resolve) {
            T::resolve(vm);
        }
    }
};

/**
 * Concrete module
 */
template <typename T>
class typed_builder
  : public typed_builder_base<T>
{
public:
    typedef typed_builder_base<T> _Base;
    typedef typename _Base::_Module_base _Base_type;

    /**
     * returns singleton instance
     */
    shared_ptr<_Base_type> fetch(po::options const& vm)
    {
        // This attaches the module-specific option values to the
        // global option values by setting an internal pointer.
        // We choose this order so subsequent module fetches from
        // the constructor of the module will get a reference to
        // the global map, and not the module-specific map.

        po::options vm_(vm);
        vm_.next(this->vm.get());

        // We use an observing weak pointer instead of an owning
        // shared pointer to let the caller decide when the
        // singleton instance and its dependencies are destroyed.
        //
        // Special care has to be taken not to destroy the
        // instance before returning it over to the caller.

        shared_ptr<T> singleton(singleton_.lock());
        if (!singleton) {
            LOG_DEBUG("instantiate module " << name());
            singleton.reset(new T(vm_));
            singleton_ = singleton;
            LOG_DEBUG("constructed module " << name());
        }
        return singleton;
    }

    /**
     * assemble module options
     */
    void options(po::options_description& desc)
    {
        _Base::options(desc);
    }

    /**
     * resolve module dependencies
     */
    void resolve(po::options const& vm)
    {
        po::options vm_(vm);
        vm_.next(this->vm.get());
        _Base::resolve(vm_);
    }

    /**
     * returns module runtime type
     */
    std::type_info const& type()
    {
        return typeid(T);
    }

    /**
     * return (demangled) module name
     */
    std::string name()
    {
        return demangled_name<T>();
    }

    /** module instance observer */
    static weak_ptr<T> singleton_;
};

template <typename T> weak_ptr<T> typed_builder<T>::singleton_;

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_BUILDER_HPP */
