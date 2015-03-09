/*
 * Copyright © 2015 Nicolas Höft
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

#ifndef TEST_UNIT_ALGORITHM_GPU_COPY_IF_KERNEL_HPP
#define TEST_UNIT_ALGORITHM_GPU_COPY_IF_KERNEL_HPP

#include <halmd/config.hpp>

template <typename T>
struct select_all
{
    HALMD_GPU_ENABLED bool operator()(T const& value) const
    {
        return true;
    }
};

template <typename T>
struct select_none
{
    HALMD_GPU_ENABLED bool operator()(T const& value) const
    {
        return false;
    }
};

template <typename T>
struct select_even
{
    HALMD_GPU_ENABLED bool operator()(T const& value) const
    {
        return (value % 2) == 0;
    }
};

template <typename T>
struct select_odd
{
    HALMD_GPU_ENABLED bool operator()(T const& value) const
    {
        return (value % 2) == 1;
    }
};

template <typename T>
struct select_prime
{
    HALMD_GPU_ENABLED bool operator()(T const& value) const
    {
        return (value % 31) == 1;   // warp size - 1
    }
};

#endif // TEST_UNIT_ALGORITHM_GPU_COPY_IF_KERNEL_HPP
