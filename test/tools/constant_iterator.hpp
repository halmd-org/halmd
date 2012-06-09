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

#ifndef TEST_TOOLS_CONSTANT_ITERATOR_HPP
#define TEST_TOOLS_CONSTANT_ITERATOR_HPP

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

/**
 * Functor that returns constant value.
 */
template <typename T>
class constant_value
{
public:
    constant_value(T const& value) : value_(value) {}

    T operator()(std::size_t) const
    {
        return value_;
    }

private:
    T value_;
};

template <typename T>
inline boost::transform_iterator<constant_value<T>, boost::counting_iterator<std::size_t>>
make_constant_iterator(T const& value, std::size_t offset)
{
    return boost::make_transform_iterator(boost::make_counting_iterator(offset), constant_value<T>(value));
}

#endif /* TEST_TOOLS_CONSTANT_ITERATOR_HPP */
