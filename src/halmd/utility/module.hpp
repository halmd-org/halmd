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

#include <boost/shared_ptr.hpp>

#include <halmd/options.hpp>

namespace halmd
{

/**
 * Module factory
 */
template <typename T>
class module
{
public:
    typedef boost::shared_ptr<T> pointer;

    /**
     * Returns instance
     */
    static pointer fetch(options const& vm)
    {
        return fetch_(vm, &T::create);
    }

    /**
     * Check if an instantiation a given type exists
     */
    static bool exists()
    {
        return exists_(&T::create);
    }

private:
    template <typename T_> friend class module;

    template <typename T_>
    static pointer fetch_(options const& vm, boost::shared_ptr<T_> (*create)(options const&))
    {
        if (!module<T_>::singleton_) {
            module<T_>::singleton_ = create(vm);
        }
        return boost::dynamic_pointer_cast<T>(module<T_>::singleton_);
    }

    template <typename T_>
    static bool exists_(boost::shared_ptr<T_> (*create)(options const&))
    {
        return module<T_>::singleton_;
    }

    static pointer singleton_;
};

template <typename T> typename module<T>::pointer module<T>::singleton_;

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_HPP */
