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

#ifndef HALMD_UTILITY_DETAIL_FACTORY_HPP
#define HALMD_UTILITY_DETAIL_FACTORY_HPP

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <exception>
#include <set>
#include <typeinfo>

#include <halmd/utility/detail/builder.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/util/logger.hpp>

namespace halmd
{
namespace utility { namespace detail
{

// import into namespace
using boost::dynamic_pointer_cast;
using boost::shared_ptr;
using boost::weak_ptr;

/**
 * A factory is implicitly instantiated once per base type
 * and holds a module singleton instance of that type.
 */
template <typename _Base>
class factory
{
public:
    typedef std::set<shared_ptr<builder<_Base> >, _builder_order<_Base> > builder_set;

    /**
     * returns module singleton instance
     */
    static shared_ptr<_Base> fetch(po::options const& vm)
    {
        // We use an observing weak pointer instead of an owning
        // shared pointer to let the caller decide when the
        // singleton instance and its dependencies are destroyed.
        //
        // Special care has to be taken not to destroy the
        // instance before returning it over to the caller.

        shared_ptr<_Base> singleton(singleton_.lock());
        if (!singleton) {
            if (builders().empty()) {
                throw std::logic_error("no modules available [" + std::string(typeid(_Base).name()) + "]");
            }
            singleton_ = singleton = (*builders().begin())->_create(vm);
        }
        return singleton;
    }

    /**
     * register module builder
     */
    static void _register(shared_ptr<builder<_Base> > builder_)
    {
        if (!builders().insert(builder_).second) {
            throw std::logic_error("module already registered [" + std::string(typeid(_Base).name()) + "]");
        }
    }

    /**
     * returns singleton builder set
     */
    static builder_set& builders()
    {
        //
        // What's the "static initialization order fiasco"?
        // http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12
        //
        static shared_ptr<builder_set> _(new builder_set);
        return *_;
    }

private:
    /** module singleton */
    static weak_ptr<_Base> singleton_;
};

template <typename _Base> weak_ptr<_Base> factory<_Base>::singleton_;

}} // namespace utility::detail

} // namespace halmd

#endif /* ! HALMD_UTILITY_DETAIL_FACTORY_HPP */
