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

#include <halmd/utility/detail/module.hpp>

namespace halmd
{

// import into top-level namespace
using utility::detail::module_exception;

/**
 * Concrete module
 */
template <typename T = void>
class module
{
public:
    typedef typename T::module_type _Base;
    typedef utility::detail::module<T> _Module;
    typedef utility::detail::builder<T> builder;
    typedef utility::detail::builder<_Base> _Base_builder;
    typedef utility::detail::factory<_Base> factory;
    typedef typename factory::builder_set builder_set;
    typedef typename builder_set::iterator builder_iterator;

    /**
     * returns singleton instance
     */
    static boost::shared_ptr<T> fetch(po::options const& vm)
    {
        LOG_DEBUG("fetch module " << typeid(T).name());
        return boost::dynamic_pointer_cast<T>(factory::fetch(vm));
    }

    static void resolve(po::options const& vm);

    /**
     * returns module name
     */
    static std::string name()
    {
        return typeid(T).name();
    }

private:
    struct _register
    {
        _register()
        {
            factory::_register(boost::shared_ptr<_Base_builder>(new _Module));
        }
    };

    static _register register_;
};

template <typename T> typename module<T>::_register module<T>::register_;

/**
 * resolve dependencies for given module
 */
template <typename T>
void module<T>::resolve(po::options const& vm)
{
    LOG_DEBUG("resolve builder " << typeid(T).name());
    builder_set& builders = factory::builders();

    for (builder_iterator it = builders.begin(); it != builders.end(); ) {
        if (!boost::dynamic_pointer_cast<builder>(*it)) {
            // module does not implement builder specification
            builders.erase(it++);
        }
        else {
            try {
                (*it)->resolve(vm);
                // resolvable builder
                return;
            }
            catch (module_exception const& e) {
                // irresolvable builder
                LOG_DEBUG(e.what());
                builders.erase(it++);
            }
        }
    }
    // no suitable modules available
    throw module_exception("irresolvable module " + module<T>::name());
}

/**
 * Type-independent module interface
 */
template <>
class module<>
{
public:
    /**
     * returns options of resolved modules
     */
    static po::options_description options()
    {
        po::options_description desc;
        utility::detail::factory<>::options(desc);
        return desc;
    }
};

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_HPP */
