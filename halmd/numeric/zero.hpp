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

#ifndef HALMD_NUMERIC_ZERO_HPP
#define HALMD_NUMERIC_ZERO_HPP

#include <halmd/config.hpp>

namespace halmd
{

/**
 * Helper to initialize arbitrary value to zero.
 *
 * This template may be specialized for types without a constructor
 * that accepts scalar values, e.g. POD (plain old datatype) structs.
 */
template <typename T>
struct zero
{
    HALMD_GPU_ENABLED operator T() const
    {
        return T(0);
    }
};

} // namespace halmd

#endif /* ! HALMD_NUMERIC_ZERO_HPP */
