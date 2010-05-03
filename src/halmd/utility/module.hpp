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

#ifndef HALMD_UTILITY_MODULE_HPP
#define HALMD_UTILITY_MODULE_HPP

#include <halmd/utility/module/builder.hpp>
#include <halmd/utility/module/exception.hpp>
#include <halmd/utility/module/factory.hpp>
#include <halmd/utility/module/module.hpp>
#include <halmd/utility/module/resolver.hpp>

namespace halmd
{

// import into top-level namespace
using utility::module::module_exception;

/**
 * Concrete module
 */
template <typename T = void>
class module
{
public:
    typedef typename T::module_type _Base;
    typedef utility::module::module<T> _Module;
    typedef utility::module::resolver<T> _Resolver;
    typedef utility::module::factory<_Base> _Base_factory;
    typedef utility::module::builder<_Base> _Base_builder;
    typedef typename _Base_factory::_Base_builder_set _Base_builder_set;
    typedef typename _Base_builder_set::iterator builder_iterator;

    /**
     * returns singleton instance
     */
    static typename _Resolver::_fetch fetch(po::options const& vm)
    {
        return _Resolver::fetch(vm);
    }

    /**
     * returns module name
     */
    static std::string name()
    {
        return _Resolver::name();
    }

    static void required(po::options const& vm);
    static void optional(po::options const& vm);
    static void many(po::options const& vm);

private:
    struct _register
    {
        _register()
        {
            _Base_factory::_register(boost::shared_ptr<_Base_builder>(new _Module));
        }
    };

    static _register register_;
};

template <typename T> typename module<T>::_register module<T>::register_;

/**
 * resolve required dependencies for given module
 */
template <typename T>
void module<T>::required(po::options const& vm)
{
    size_t count = _Resolver::resolve(vm);

    if (count == 0) {
        throw module_exception("irresolvable dependency " + module<T>::name());
    }
    if (count > 1) {
        throw module_exception("ambiguous dependency " + module<T>::name());
    }
}

/**
 * resolve optional dependencies for given module
 */
template <typename T>
void module<T>::optional(po::options const& vm)
{
    size_t count = _Resolver::resolve(vm);

    if (count > 1) {
        throw module_exception("ambiguous dependency " + module<T>::name());
    }
}

/**
 * resolve one-to-many dependencies for given module
 */
template <typename T>
void module<T>::many(po::options const& vm)
{
    _Resolver::resolve(vm);
}

/**
 * Type-independent module interface
 */
template <>
class module<>
{
public:
    typedef utility::module::resolver<> _Resolver;

    /**
     * returns options of resolved modules
     */
    static po::options_description options()
    {
        return _Resolver::options();
    }
};

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_HPP */
