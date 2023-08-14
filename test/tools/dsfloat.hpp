/*
 * Copyright © 2017 Daniel Kirchner
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

#ifndef TEST_TOOLS_DSFLOAT_HPP
#define TEST_TOOLS_DSFLOAT_HPP

#include <halmd/numeric/mp/dsfloat.hpp>
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <limits>
#include <ostream>

namespace boost {
namespace math {
namespace fpc {

template <>
struct tolerance_based<halmd::dsfloat> : boost::true_type{};

} // namespace fpc
} // namespace math
} // namespace boost

// FIXME define numeric_limits for dsfloat
// see, e.g., http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
template<typename T>
struct dsfloat_aware_numeric_limits : public std::numeric_limits<T> {
};
template<>
struct dsfloat_aware_numeric_limits<halmd::dsfloat> {
    static constexpr double epsilon() { return 5.6843418860808015e-14;  } // std::pow(double(2), -44);
    static constexpr float min() { return std::numeric_limits<float>::min(); }
};

#endif // TEST_TOOLS_DSFLOAT_HPP
