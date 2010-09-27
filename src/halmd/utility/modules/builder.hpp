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

#ifndef HALMD_UTILITY_MODULES_BUILDER_HPP
#define HALMD_UTILITY_MODULES_BUILDER_HPP

#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_object.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/weak_ptr.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/options.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/modules/concept.hpp>
#include <halmd/utility/modules/exception.hpp>

namespace halmd
{
namespace modules
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
    virtual void depends() = 0;
    virtual void select(po::variables_map const& vm) = 0;
    virtual std::string name() = 0;
    virtual std::type_info const& type() = 0;
    po::variables_map vm;
};

/**
 * A base template, which corresponds to the base class of a
 * module (e.g. force, particle). The module fetch function
 * casts type-independent builder pointers to this type to
 * create module instances.
 */
template <typename T, typename Factory, typename Enable = void>
struct typed_builder_base
  : public untyped_builder_base
{
    typedef T BaseT;
    virtual shared_ptr<BaseT> fetch(Factory& factory, po::variables_map const& vm) = 0;

    /**
     * resolve class dependencies
     */
    virtual void depends()
    {
        T::depends();
    }

    /**
     * resolve class dependencies
     */
    virtual void select(po::variables_map const& vm)
    {
        T::select(vm);
    }
};

/**
 * A deriving template to the base template.
 */
template <typename T, typename Factory>
struct typed_builder_base<T, Factory, typename enable_if<is_object<typename T::_Base> >::type>
  : public typed_builder_base<typename T::_Base, Factory>
{
    typedef typed_builder_base<typename T::_Base, Factory> Base;
    typedef typename Base::BaseT BaseT;

    /**
     * resolve class dependencies
     */
    virtual void depends()
    {
        Base::depends();
        // check for inherited base module function
        boost::function_requires<modules::ModuleConcept<T> >();
        T::depends();
    }

    /**
     * resolve class dependencies
     */
    virtual void select(po::variables_map const& vm)
    {
        Base::select(vm);
        // check for inherited base module function
        boost::function_requires<modules::ModuleConcept<T> >();
        T::select(vm);
    }
};

/**
 * Concrete module
 */
template <typename T, typename Factory>
class typed_builder
  : public typed_builder_base<T, Factory>
{
public:
    typedef typed_builder_base<T, Factory> Base;
    typedef typename Base::BaseT BaseT;

    /**
     * returns singleton instance
     */
    shared_ptr<BaseT> fetch(Factory& factory, po::variables_map const& vm)
    {
        // This attaches the module-specific option values to the
        // global option values by setting an internal pointer.
        // We choose this order so subsequent module fetches from
        // the constructor of the module will get a reference to
        // the global map, and not the module-specific map.

        po::variables_map vm_(vm);
        vm_.next(&this->vm);

        // We use an observing weak pointer instead of an owning
        // shared pointer to let the caller decide when the
        // singleton instance and its dependencies are destroyed.
        //
        // Special care has to be taken not to destroy the
        // instance before returning it over to the caller.

        shared_ptr<T> singleton(singleton_.lock());
        if (!singleton) {
            LOG_DEBUG("instantiate module " << name());
            singleton.reset(new T(factory, vm_));
            singleton_ = singleton;
            LOG_DEBUG("constructed module " << name());
        }
        return singleton;
    }

    /**
     * resolve module dependencies
     */
    void select(po::variables_map const& vm)
    {
        po::variables_map vm_(vm);
        vm_.next(&this->vm);
        Base::select(vm_);
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
    weak_ptr<T> singleton_;
};

} // namespace modules

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULES_BUILDER_HPP */
