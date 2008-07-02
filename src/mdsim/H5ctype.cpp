/* C types to HDF5 native data type translation
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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
#include "H5ctype.hpp"

template <> H5::PredType const& H5ctype<int8_t>::type = H5::PredType::NATIVE_INT8;
template <> H5::PredType const& H5ctype<uint8_t>::type = H5::PredType::NATIVE_UINT8;
template <> H5::PredType const& H5ctype<int16_t>::type = H5::PredType::NATIVE_INT16;
template <> H5::PredType const& H5ctype<uint16_t>::type = H5::PredType::NATIVE_UINT16;
template <> H5::PredType const& H5ctype<int32_t>::type = H5::PredType::NATIVE_INT32;
template <> H5::PredType const& H5ctype<uint32_t>::type = H5::PredType::NATIVE_UINT32;
template <> H5::PredType const& H5ctype<int64_t>::type = H5::PredType::NATIVE_INT64;
template <> H5::PredType const& H5ctype<uint64_t>::type = H5::PredType::NATIVE_UINT64;
template <> H5::PredType const& H5ctype<float>::type = H5::PredType::NATIVE_FLOAT;
template <> H5::PredType const& H5ctype<double>::type = H5::PredType::NATIVE_DOUBLE;
