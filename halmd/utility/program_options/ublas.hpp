/*
 * Copyright © 2011  Peter Colberg
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

#ifndef HALMD_UTILITY_PROGRAM_OPTIONS_UBLAS_HPP
#define HALMD_UTILITY_PROGRAM_OPTIONS_UBLAS_HPP

#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>

/**
 * This file adds support for parsing boost::numeric::ublas::matrix
 * and boost::numeric::ublas::vector as option values using Boost
 * program_options.
 *
 * A vector is represented as a comma-delimited string, e.g. 1,2,3,4
 *
 * A matrix is represented as colon-delimited rows with each row as a
 * comma-delimited string, e.g. 11,12,13:21,22,23:31,32,33
 */
namespace std {

/**
 * Read Boost uBLAS matrix from input stream
 */
template <typename T>
std::istream& operator>>(std::istream& is, boost::numeric::ublas::matrix<T>& value)
{
    namespace ublas = boost::numeric::ublas;
    std::vector<ublas::vector<T> > row;
    while (!is.eof()) {
        std::string str;
        getline(is, str, ':');
        row.push_back(boost::lexical_cast<ublas::vector<T> >(str));
    }
    if (!row.empty()) {
        value.resize(row.size(), row.front().size());
        for (std::size_t i = 0; i < row.size(); ++i) {
            if (!(row[i].size() == value.size2())) {
                throw boost::bad_lexical_cast();
            }
            ublas::matrix_row<ublas::matrix<T> >(value, i) = row[i];
        }
    }
    return is;
}

/**
 * Write Boost uBLAS matrix to output stream
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, boost::numeric::ublas::matrix<T> const& value)
{
    namespace ublas = boost::numeric::ublas;
    for (std::size_t i = 0; i < value.size1(); ++i) {
        if (i > 0) {
            os << ':';
        }
        ublas::vector<T> const& row = ublas::matrix_row<ublas::matrix<T> const>(value, i);
        os << row;
    }
    return os;
}

/**
 * Read Boost uBLAS vector from input stream
 */
template <typename T>
std::istream& operator>>(std::istream& is, boost::numeric::ublas::vector<T>& value)
{
    std::vector<T> v;
    while (!is.eof()) {
        std::string str;
        getline(is, str, ',');
        v.push_back(boost::lexical_cast<T>(str));
    }
    value.resize(v.size());
    copy(v.begin(), v.end(), value.begin());
    return is;
}

/**
 * Write Boost uBLAS vector to output stream
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, boost::numeric::ublas::vector<T> const& value)
{
    for (std::size_t i = 0; i < value.size(); ++i) {
        if (i > 0) {
            os << ',';
        }
        os << value(i);
    }
    return os;
}

} // namespace std

#endif /* ! HALMD_UTILITY_PROGRAM_OPTIONS_UBLAS_HPP */
