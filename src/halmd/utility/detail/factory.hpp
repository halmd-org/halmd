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
#include <exception>
#include <set>
#include <typeinfo>

#include <halmd/utility/detail/builder.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace utility { namespace detail
{

template <typename T>
struct factory
{
    typedef boost::shared_ptr<T> base_ptr;
    typedef boost::shared_ptr<detail::builder<T> > builder_base_ptr;
    typedef typename detail::builder<T>::less builder_compare;
    typedef std::set<builder_base_ptr, builder_compare> builder_set;

    static base_ptr fetch(po::options const& vm)
    {
        if (!singleton) {
            singleton = builder->create(vm);
        }
        return singleton;
    }

    static void register_(builder_base_ptr const& builder)
    {
        std::string const type = typeid(T).name();
        if (!builders()->insert(builder).second) {
            throw std::logic_error("module already registered [" + type + "]");
        }
    }

    //
    // What's the "static initialization order fiasco"?
    // http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12
    //
    static boost::shared_ptr<builder_set> builders()
    {
        static boost::shared_ptr<builder_set> _(new builder_set);
        return _;
    }

    static base_ptr singleton;
    static builder_base_ptr builder;
};

template <typename T> typename factory<T>::base_ptr factory<T>::singleton;
template <typename T> typename factory<T>::builder_base_ptr factory<T>::builder;

}} // namespace utility::detail

} // namespace halmd

#endif /* ! HALMD_UTILITY_DETAIL_FACTORY_HPP */
