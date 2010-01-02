/* HDF5 C++ extensions
 *
 * Copyright © 2008-2009  Peter Colberg
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

#include <stdint.h>

#include <halmd/util/H5xx.hpp>

namespace H5xx
{

template <> H5::PredType const& ctype<int8_t>::type = H5::PredType::NATIVE_INT8;
template <> H5::PredType const& ctype<uint8_t>::type = H5::PredType::NATIVE_UINT8;
template <> H5::PredType const& ctype<int16_t>::type = H5::PredType::NATIVE_INT16;
template <> H5::PredType const& ctype<uint16_t>::type = H5::PredType::NATIVE_UINT16;
template <> H5::PredType const& ctype<int32_t>::type = H5::PredType::NATIVE_INT32;
template <> H5::PredType const& ctype<uint32_t>::type = H5::PredType::NATIVE_UINT32;
template <> H5::PredType const& ctype<int64_t>::type = H5::PredType::NATIVE_INT64;
template <> H5::PredType const& ctype<uint64_t>::type = H5::PredType::NATIVE_UINT64;
template <> H5::PredType const& ctype<float>::type = H5::PredType::NATIVE_FLOAT;
template <> H5::PredType const& ctype<double>::type = H5::PredType::NATIVE_DOUBLE;

} // namespace H5xx
