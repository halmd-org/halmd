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

#ifndef HALMD_UTILITY_PROGRAM_OPTIONS_ARRAY_HPP
#define HALMD_UTILITY_PROGRAM_OPTIONS_ARRAY_HPP

#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <vector>

namespace std {

/**
 * extract comma-separated option values into fixed-size array
 */
template <typename T, size_t size>
std::istream& operator>>(std::istream& is, boost::array<T, size>& value)
{
    BOOST_FOREACH(T& v, value) {
        std::string str;
        getline(is, str, ',');
        v = boost::lexical_cast<T>(str);
    }
    return is;
}

template <typename T, size_t size>
std::ostream& operator<<(std::ostream& os, boost::array<T, size> const& value)
{
    BOOST_FOREACH(T const& v, value) {
        if (&v != &value.front()) {
            os << ',';
        }
        os << v;
    }
    return os;
}

/**
 * extract comma-separated option values into variable-size array
 */
template <typename T>
std::istream& operator>>(std::istream& is, boost::multi_array<T, 1>& value)
{
    std::vector<T> v;
    std::string str;
    while (!is.eof()) {
        getline(is, str, ',');
        v.push_back(boost::lexical_cast<T>(str));
    }
    boost::array<size_t, 1> extents = {{ v.size() }};
    value.resize(extents);
    std::copy(v.begin(), v.end(), value.begin());
    return is;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, boost::multi_array<T, 1> const& value)
{
    BOOST_FOREACH(T const& v, value) {
        if (&v != &(*value.begin())) {
            os << ',';
        }
        os << v;
    }
    return os;
}

} // namespace std

#endif /* ! HALMD_UTILITY_PROGRAM_OPTIONS_ARRAY_HPP */
