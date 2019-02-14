/*
 * Copyright © 2013 Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_UTILITY_OWNER_EQUAL_HPP
#define HALMD_UTILITY_OWNER_EQUAL_HPP

#include <memory>

namespace halmd {

/**
 * Provide owner_equal(a,b) for std::shared_ptr and std::weak_ptr. It returns
 * true if the two arguments refer to the same managed object.
 *
 * Based on http://stackoverflow.com/questions/12301916/equality-compare-stdweak-ptr.
 */

template <typename T, typename U>
inline bool owner_equal(std::weak_ptr<T> const& t, std::weak_ptr<U> const& u)
{
    return !t.owner_before(u) && !u.owner_before(t);
}

template <typename T, typename U>
inline bool owner_equal(std::shared_ptr<T> const& t, std::shared_ptr<U> const& u)
{
    return !t.owner_before(u) && !u.owner_before(t);
}

template <typename T, typename U>
inline bool owner_equal(std::weak_ptr<T> const& t, std::shared_ptr<U> const& u)
{
    return !t.owner_before(u) && !u.owner_before(t);
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_OWNER_EQUAL_HPP */
