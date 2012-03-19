/*
 * Copyright © 2012  Peter Colberg
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

#ifndef HALMD_NUMERIC_CAST_HPP
#define HALMD_NUMERIC_CAST_HPP

#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <stdexcept>

namespace halmd {

/**
 * This function validates at runtime whether the value is exactly
 * representable in type T, and otherwise throws an exception.
 */
template <typename T, typename S>
inline typename boost::enable_if<boost::is_convertible<S, T>, T>::type
checked_narrowing_cast(S const& value)
{
    T result = T(value);
    if (S(result) != value) {
        throw std::domain_error("value not exactly representable in target type");
    }
    return result;
}

} // namespace halmd

#endif /* ! HALMD_NUMERIC_CAST_HPP */
