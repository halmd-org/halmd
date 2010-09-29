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

#ifndef HALMD_UTILITY_READ_INTEGER_HPP
#define HALMD_UTILITY_READ_INTEGER_HPP

#include <boost/type_traits/is_integral.hpp>
#include <boost/utility/enable_if.hpp>
#include <fstream>
#include <stdexcept>

namespace halmd
{

/**
 * Read number from file
 */
template <typename T>
typename boost::enable_if<boost::is_integral<T>, T>::type
read_integer(std::string const& file_name)
{
    T value;
    try {
        std::ifstream f;
        f.exceptions(std::ifstream::eofbit|std::ifstream::failbit|std::ifstream::badbit);
        f.open(file_name.c_str());
        f.read(reinterpret_cast<char*>(&value), sizeof(value));
        f.close();
    }
    catch (std::ifstream::failure const& e) {
        throw std::runtime_error("failed to read integer from file: " + file_name);
    }
    return value;
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_READ_INTEGER_HPP */
