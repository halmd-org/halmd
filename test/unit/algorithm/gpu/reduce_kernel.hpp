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

#ifndef TEST_UNIT_ALGORITHM_GPU_REDUCE_KERNEL_HPP
#define TEST_UNIT_ALGORITHM_GPU_REDUCE_KERNEL_HPP

#include <halmd/config.hpp>

template <typename input_type, typename output_type>
class sum
{
public:
    typedef input_type argument_type;

    sum(output_type const& init = 0) : sum_(init) {}

    output_type operator()() const
    {
        return sum_;
    }

    HALMD_GPU_ENABLED void operator()(argument_type const& value)
    {
        sum_ += value;
    }

    HALMD_GPU_ENABLED void operator()(sum const& acc)
    {
        sum_ += acc.sum_;
    }

private:
    output_type sum_;
};

template <typename input_type, typename output_type>
class sum_of_squares
{
public:
    typedef input_type first_argument_type;
    typedef input_type second_argument_type;

    sum_of_squares(output_type const& init = 0) : sum_(init) {}

    output_type operator()() const
    {
        return sum_;
    }

    HALMD_GPU_ENABLED void operator()(first_argument_type const& first, second_argument_type const& second)
    {
        sum_ += output_type(first) * output_type(second);
    }

    HALMD_GPU_ENABLED void operator()(sum_of_squares const& acc)
    {
        sum_ += acc.sum_;
    }

private:
    output_type sum_;
};

#endif /* ! TEST_UNIT_ALGORITHM_GPU_REDUCE_KERNEL_HPP */
