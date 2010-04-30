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
#include <halmd/utility/detail/exception.hpp>
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

template <typename _Base = void>
class factory;

/**
 * Factory registry
 */
template <>
class factory<>
{
public:
    typedef std::set<shared_ptr<factory<> > > factory_set;
    typedef factory_set::iterator factory_iterator;
    virtual ~factory() {}

    /**
     * assemble options of resolved builders
     */
    static void options(po::options_description& desc)
    {
        factory_set& factories = factory<>::factories();

        for (factory_iterator it = factories.begin(); it != factories.end(); ++it) {
            (*it)->_options(desc);
        }
    }

protected:
    virtual void _options(po::options_description& desc) = 0;

    /**
     * register factory
     */
    struct _register
    {
        _register(shared_ptr<factory<> > factory_)
        {
            factory<>::factories().insert(factory_);
        }
    };

private:
    /**
     * returns singleton factory set
     */
    static factory_set& factories()
    {
        static shared_ptr<factory_set> _(new factory_set);
        return *_;
    }
};

/**
 * A factory is implicitly instantiated once per base type
 * and holds a module singleton instance of that type.
 */
template <typename _Base>
class factory
  : public factory<>
{
public:
    typedef std::set<shared_ptr<builder<_Base> >, _builder_order<_Base> > builder_set;

    /**
     * returns module singleton instance
     */
    static shared_ptr<_Base> fetch(po::options const& vm)
    {
        if (_builders().empty()) {
            throw module_exception(std::string("unavailable module ") + typeid(_Base).name());
        }
        return (*_builders().begin())->fetch(vm);
    }

    /**
     * register module builder
     */
    static void _register(shared_ptr<builder<_Base> > builder_)
    {
        if (!_builders().insert(builder_).second) {
            throw module_exception(std::string("duplicate builder ") + typeid(_Base).name());
        }
    }

    /**
     * returns singleton builder set
     */
    static builder_set& builders()
    {
        // Calling this function signals that the factory is in use.
        // Therefore we register the factory (once) to enable options
        // collection of the builder with the highest rank.

        static factory<>::_register register_(shared_ptr<factory<> >(new factory<_Base>));
        return _builders();
    }

protected:
    /**
     * returns module options of builder with highest rank
     */
    void _options(po::options_description& desc)
    {
        if (_builders().empty()) {
            throw module_exception(std::string("unavailable module ") + typeid(_Base).name());
        }
        (*_builders().begin())->options(desc);
    }

private:
    static builder_set& _builders()
    {
        // What's the "static initialization order fiasco"?
        // http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12

        static shared_ptr<builder_set> _(new builder_set);
        return *_;
    }
};

}} // namespace utility::detail

} // namespace halmd

#endif /* ! HALMD_UTILITY_DETAIL_FACTORY_HPP */
