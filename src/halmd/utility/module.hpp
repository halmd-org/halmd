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

#include <halmd/utility/module/module.hpp>

namespace halmd
{

// import into top-level namespace
using boost::shared_ptr;
using utility::module::module;

// forward compatibility with upcoming module mechanism rewrite
namespace modules
{

template <typename T, typename U>
void required()
{
    module<U>::required(*utility::module::factory::vm);
}

template <typename T, typename U>
void optional()
{
    module<U>::optional(*utility::module::factory::vm);
}

template <typename T>
typename module<T>::_fetch fetch(po::options const& vm)
{
    return module<T>::fetch(vm);
}

} // namespace modules

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_HPP */
