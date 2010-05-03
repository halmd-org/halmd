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

#ifndef HALMD_UTILITY_MODULE_FACTORY_HPP
#define HALMD_UTILITY_MODULE_FACTORY_HPP

#include <boost/shared_ptr.hpp>
#include <set>
#include <typeinfo>

#include <halmd/utility/module/builder.hpp>
#include <halmd/utility/module/exception.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/util/logger.hpp>

namespace halmd
{
namespace utility { namespace module
{

// import into namespace
using boost::shared_ptr;

/**
 * A factory is implicitly instantiated once per base type.
 */
template <typename _Base>
class factory
{
public:
    typedef shared_ptr<builder<_Base> > _Base_builder_ptr;
    typedef std::set<_Base_builder_ptr, _builder_rank> _Base_builder_set;

    /**
     * register module builder
     */
    static void _register(_Base_builder_ptr builder_)
    {
        if (!modules().insert(builder_).second) {
            throw module_exception("duplicate builder " + name());
        }
    }

    /**
     * returns singleton builder set
     */
    static _Base_builder_set& modules()
    {
        // What's the "static initialization order fiasco"?
        // http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12

        static _Base_builder_set modules_;
        return modules_;
    }

    /**
     * returns module name
     */
    static std::string name()
    {
        return typeid(_Base).name();
    }
};

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_FACTORY_HPP */
