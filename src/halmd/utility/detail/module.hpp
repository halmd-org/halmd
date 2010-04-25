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

#ifndef HALMD_UTILITY_DETAIL_MODULE_HPP
#define HALMD_UTILITY_DETAIL_MODULE_HPP

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <exception>
#include <set>
#include <typeinfo>

#include <halmd/options.hpp>
#include <halmd/utility/detail/factory.hpp>
#include <halmd/utility/detail/builder.hpp>

namespace halmd
{
namespace utility { namespace detail
{

template <typename T>
class module
{
public:
    typedef typename T::module_type base_type;
    typedef detail::factory<base_type> factory;
    typedef boost::shared_ptr<detail::builder<base_type> > builder_base_ptr;
    typedef detail::builder<T> builder_type;

    static boost::shared_ptr<T> fetch(options const& vm)
    {
        return boost::dynamic_pointer_cast<T>(factory::fetch(vm));
    }

    static void resolve(options const& vm);

private:
    /**
     * Register module builder singleton
     */
    struct register_
    {
        register_()
        {
            factory::register_(builder_base_ptr(new builder_type));
        }
    };

    static register_ register__;
};

template <typename T> typename module<T>::register_ module<T>::register__;

/**
 * Resolve module dependency
 */
template <typename T>
void module<T>::resolve(options const& vm)
{
    if (factory::builder) {
        if (boost::dynamic_pointer_cast<builder_type>(factory::builder)) {
            return;
        }
    }
    BOOST_FOREACH(builder_base_ptr const& builder, *factory::builders()) {
        if (!boost::dynamic_pointer_cast<builder_type>(builder)) {
            continue;
        }
        if (!factory::builder || factory::builder->is_base_of(builder)) {
            try {
                builder->resolve(vm);
            }
            catch (std::exception const&) {
                continue;
            }
            factory::builder.reset(new builder_type(*builder));
            return;
        }
    }
    std::string type = typeid(T).name();
    throw std::logic_error("no modules available [" + type + "]");
}

}} // namespace utility::detail

} // namespace halmd

#endif /* ! HALMD_UTILITY_DETAIL_MODULE_HPP */
