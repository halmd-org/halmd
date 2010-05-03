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

#ifndef HALMD_UTILITY_MODULE_RESOLVER_HPP
#define HALMD_UTILITY_MODULE_RESOLVER_HPP

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <deque>
#include <set>
#include <typeinfo>
#include <vector>

#include <halmd/utility/module/builder.hpp>
#include <halmd/utility/module/factory.hpp>
#include <halmd/utility/module/exception.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/util/logger.hpp>

namespace halmd
{
namespace utility { namespace module
{

// import into namespace
using boost::dynamic_pointer_cast;
using boost::shared_ptr;

template <typename T = void>
class resolver;

/**
 * Module dependency resolution
 */
template <typename T>
class resolver
{
public:
    typedef typename T::module_type _Base;
    typedef shared_ptr<builder<_Base> > _Base_builder_ptr;
    typedef std::set<_Base_builder_ptr, _builder_rank_exclude_base> _Module_set;
    typedef factory<_Base> _Base_factory;
    typedef builder<T> _Builder;

    static size_t resolve(po::options const& vm);

    /**
     * returns module singleton instance(s)
     */
    struct _fetch
    {
        _fetch(po::options const& vm) : vm(vm) {}
        options const& vm;

        /**
         * returns required or optional instance
         */
        operator shared_ptr<T>()
        {
            shared_ptr<T> result;
            if (!modules_->empty()) {
                result = dynamic_pointer_cast<T>((*modules_->begin())->fetch(vm));
            }
            return result;
        }

        /**
         * returns many instances
         */
        operator std::vector<shared_ptr<T> >()
        {
            std::vector<shared_ptr<T> > result;
            BOOST_FOREACH(_Base_builder_ptr builder, *modules_) {
                result.push_back(dynamic_pointer_cast<T>(builder->fetch(vm)));
            }
            return result;
        }
    };

    static _fetch fetch(po::options const& vm)
    {
        if (!modules_) {
            throw module_exception("missing dependencies for module " + name());
        }
        return _fetch(vm);
    }

    /**
     * returns module name
     */
    static std::string name()
    {
        return typeid(T).name();
    }

private:
    static shared_ptr<_Module_set> modules_;
};

template <typename T> shared_ptr<typename resolver<T>::_Module_set> resolver<T>::modules_;

/**
 * A type-independent class to keep track of used modules.
 */
template <>
class resolver<>
{
public:
    typedef builder<> _Builder;
    typedef shared_ptr<_Builder > _Builder_ptr;
    typedef std::deque<_Builder_ptr> _Module_stack;

    /**
     * assemble module options
     */
    static po::options_description options()
    {
        // Using the stack of modules created during dependency
        // resolution we build a unique set of used modules and
        // collect the options of each module.

        std::set<_Builder_ptr> set(stack().begin(), stack().end());
        po::options_description desc;
        std::for_each(set.begin(), set.end(), boost::bind(&_Builder::options, _1, boost::ref(desc)));
        return desc;
    }

private:
    template <typename T> friend class resolver;

    static _Module_stack& stack()
    {
        // A static class (not template) member variable has to
        // be defined in a separate source file, which we avoid
        // by using a static local variable instead.

        static _Module_stack stack_;
        return stack_;
    }
};

/**
 * resolve module dependencies
 *
 * returns the number of suitable modules, to can be used to
 * validate a dependency as required, optional or one-to-many.
 */
template <typename T>
size_t resolver<T>::resolve(po::options const& vm)
{
    // Given a required, optional or one-to-many dependency,
    // which may be a (multiply) derived or base class, this
    // function looks for suitable modules.
    // A module is suitable (1) if it is or derives from the
    // given dependency type and (2) if its dependencies can be
    // resolved.
    //

    resolver<>::_Module_stack& stack = resolver<>::stack();
    typedef resolver<>::_Module_stack::iterator stack_iterator;

#define _LOG_DEBUG_STACK(x) LOG_DEBUG("[" << stack.size() << "] " << x)

    // We cache the result of a dependency resolution by using a
    // shared pointer to the set of suitable modules and checking
    // the pointer for validity.

    if (!modules_) {
        modules_.reset(new _Module_set);

        _LOG_DEBUG_STACK("resolve dependency " << name());

        // Check each of the modules registered in the base class
        // factory during program startup for suitability.

        BOOST_FOREACH(_Base_builder_ptr module, _Base_factory::modules()) {
            // Check if the module is or derives from the given
            // dependency type with a downcast to this type.
            if (!dynamic_pointer_cast<_Builder>(module)) {
                continue;
            }
            // Check if the set of suitable modules contains a
            // module which derives from this module. The module
            // set ordering guarantees that derived modules
            // come before base modules in the sequence.
            if (modules_->count(module)) {
                // Derived modules should override base modules,
                // therefore we skip this base module.
                continue;
            }

            _LOG_DEBUG_STACK("resolve module " << module->name());

            // Take note of the current top of the dependency
            // resolution stack to rewind in case of failure. 
            stack_iterator top = stack.end();
            stack.push_back(module);
            // Try to resolve the dependencies of the module.
            try {
                module->resolve(vm);
            }
            catch (module_exception const& e) {
                // The module is irresolvable, therefore we
                // rewind the stack and skip this module.
                _LOG_DEBUG_STACK(e.what());
                stack.erase(top, stack.end());
                continue;
            }
            // The module is resolvable.
            modules_->insert(module);
        }
    }

#undef _LOG_DEBUG_STACK

    return modules_->size();
}

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_RESOLVER_HPP */
