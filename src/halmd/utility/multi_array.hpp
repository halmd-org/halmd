/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_UTILITY_MULTI_ARRAY_HPP
#define HALMD_UTILITY_MULTI_ARRAY_HPP

#include <algorithm>
#include <boost/multi_array.hpp>
#include <vector>

namespace halmd
{

/**
 * Make 1-dimensional boost::multi_array from std::vector
 */
template <typename T>
boost::multi_array<T, 1> make_multi_array(std::vector<T> const& vector)
{
    boost::multi_array<T, 1> array(boost::extents[vector.size()]);
    std::copy(vector.begin(), vector.end(), array.begin());
    return array;
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_MULTI_ARRAY_HPP */
