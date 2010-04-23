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
#include <exception>
#include <list>
#include <typeinfo>

#include <halmd/options.hpp>

namespace halmd
{

// forward declaration
template <typename T>
class module;

/**
 * Module factory
 */
template <typename T>
class factory
{
public:
    typedef boost::shared_ptr<T> pointer;
    typedef boost::shared_ptr<T> (*function)(options const&);

private:
    template <typename T_> friend class module;

    typedef std::list<function> functions;
    typedef boost::shared_ptr<functions> functions_ptr;

    //
    // What's the "static initialization order fiasco"?
    // http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12
    //
    static functions_ptr modules()
    {
        static functions_ptr modules_(new functions);
        return modules_;
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
        BOOST_FOREACH (function const& create, *modules()) {
            singleton_ = create(vm);
            if (singleton_) {
                return singleton_;
            }
        }
        std::string type = typeid(T).name();
        throw std::logic_error("no modules registered [" + type + "]");
    }
    return singleton_;
}

/**
 * Module
 */
template <typename T>
class module
{
public:
    typedef boost::shared_ptr<T> pointer;

    /**
     * Returns singleton instance of derived type
     */
    static pointer fetch(options const& vm)
    {
        return boost::dynamic_pointer_cast<T>(factory_::fetch(vm));
    }

private:
    typedef typename T::factory_ factory_;

    /**
     * Register module singleton
     */
    struct register_
    {
        register_()
        {
            factory_::modules()->push_back(&T::create);
        }
    };

    static register_ register__;
};

template <typename T> typename module<T>::register_ module<T>::register__;

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_HPP */
