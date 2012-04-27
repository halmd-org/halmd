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

#include <cstddef> // std::size_t

/**
 * Iterator for BOOST_CHECK_EQUAL_COLLECTIONS that returns constant.
 *
 * Note that this is *not* a complete iterator with regard to any of the
 * iterator categories of the STL, merely a helper for comparing a range
 * of values with a constant using BOOST_CHECK_EQUAL_COLLECTIONS.
 */
template <typename T>
class constant_iterator
{
public:
    constant_iterator(T const& value, std::size_t index) : value_(value), index_(index) {}

    bool operator!=(constant_iterator const& other)
    {
        return index_ != other.index_;
    }

    constant_iterator& operator++()
    {
        ++index_;
        return *this;
    }

    T const& operator*() const
    {
        return value_;
    }

private:
    T value_;
    std::size_t index_;
};

#endif /* TEST_TOOLS_CONSTANT_ITERATOR_HPP */
