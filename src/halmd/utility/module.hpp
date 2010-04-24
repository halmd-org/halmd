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

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <exception>
#include <set>
#include <typeinfo>

#include <halmd/options.hpp>

namespace halmd
{

// forward declaration
template <typename T>
class module;

/**
 * FIXME
 */
template <typename T>
class factory
{
public:
    typedef boost::shared_ptr<T> pointer;

private:
    template <typename T_> friend class module;

    typedef std::set<module<T> > module_set;

    //
    // What's the "static initialization order fiasco"?
    // http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12
    //
    static boost::shared_ptr<module_set> modules()
    {
        static boost::shared_ptr<module_set> modules_(new module_set);
        return modules_;
    }

    static void register_(module<T> const& module_)
    {
        std::string const type = typeid(T).name();
        if (!modules()->insert(module_).second) {
            throw std::logic_error("module already registered [" + type + "]");
        }
    }

    static pointer fetch(options const& vm);
    static pointer singleton_;
};

template <typename T> typename factory<T>::pointer factory<T>::singleton_;

/**
 * Returns singleton instance of base type
 */
template <typename T>
typename factory<T>::pointer factory<T>::fetch(options const& vm)
{
    if (!singleton_) {
        BOOST_FOREACH (module<T> const& module_, *modules()) {
            singleton_ = module_.create_(vm);
            if (singleton_) {
                return singleton_;
            }
        }
        std::string type = typeid(T).name();
        throw std::logic_error("no matching modules available [" + type + "]");
    }
    return singleton_;
}

/**
 * FIXME
 */
template <typename T, typename Enable = void>
struct _module_priority
{
    enum { value = _module_priority<typename T::_Base>::value + 1 };
};

template <typename T>
struct _module_priority<T, typename boost::enable_if<boost::is_same<T, typename T::module_ptr::value_type> >::type>
{
    enum { value = 0 };
};

/**
 * FIXME
 */
template <typename T>
class module
{
public:
    typedef boost::shared_ptr<T> pointer;

    module(boost::shared_ptr<T> (*create)(options const&), int priority)
      : create_(create)
      , priority_(priority)
    {}

    bool operator<(module<T> const& module_) const
    {
        if (priority_ == module_.priority_) {
            // order of modules with same priority is undefined
            return create_ != module_.create_;
        }
        // order module with higher priority *before* module with lower priority
        return priority_ > module_.priority_;
    }

    /**
     * Returns singleton instance of derived type
     */
    static pointer fetch(options const& vm)
    {
        typedef typename T::module_ptr::value_type _Base;
        return boost::dynamic_pointer_cast<T>(factory<_Base>::fetch(vm));
    }

private:
    template <typename T_> friend class factory;

    /**
     * Register module singleton
     */
    struct register_
    {
        register_()
        {
            register__(&T::create);
        }

        template <typename T_>
        void register__(boost::shared_ptr<T_> (*create)(options const&)) {
            factory<T_>::register_(module<T_>(&T::create, _module_priority<T>::value));
        }
    };

    static register_ register__;
    boost::shared_ptr<T> (*create_)(options const&);
    int priority_;
};

template <typename T> typename module<T>::register_ module<T>::register__;

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_HPP */
