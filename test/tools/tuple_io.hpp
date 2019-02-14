/*
 * Copyright © 2014 Nicolas Höft
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

#ifndef HALMD_TEST_TOOLS_TUPLE_IO_HPP
#define HALMD_TEST_TOOLS_TUPLE_IO_HPP

#include <ostream>
#include <tuple>

namespace tuple_io {

// Define a type which holds an unsigned integer value
template<std::size_t> struct int_{};

template <class Tuple, size_t Pos>
std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<Pos> ) {
  out << std::get<std::tuple_size<Tuple>::value - Pos>(t) << ' ';
  return print_tuple(out, t, int_<Pos - 1>());
}

template <class Tuple>
std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<1> ) {
  return out << std::get<std::tuple_size<Tuple>::value - 1>(t);
}

} // namespace tuple_io

template <class... Args>
std::ostream& operator<<(std::ostream& out, const std::tuple<Args...>& t) {
  out << '(';
  tuple_io::print_tuple(out, t, tuple_io::int_<sizeof...(Args)>());
  return out << ')';
}

#endif /* ! HALMD_TEST_TOOLS_TUPLE_IO_HPP */
