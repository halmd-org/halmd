/*
 * Copyright © 2012 Peter Colberg
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

#ifndef TEST_UNIT_ALGORITHM_GPU_REDUCE_KERNEL_HPP
#define TEST_UNIT_ALGORITHM_GPU_REDUCE_KERNEL_HPP

#include <halmd/config.hpp>
#include <halmd/utility/iterator.hpp>

template <typename input_type, typename output_type>
class sum
{
public:
    typedef input_type const* iterator;

    // DefaultConstructible
    explicit sum(output_type const& init = 0) : sum_(init) {}

    output_type operator()() const
    {
        return sum_;
    }

    HALMD_GPU_ENABLED void operator()(input_type const& value)
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
class sum_with_constant
{
public:
    typedef halmd::zip_iterator<input_type const*, halmd::constant_iterator<input_type> > iterator;

    // DefaultConstructible
    explicit sum_with_constant(output_type const& init = 0) : sum_(init) {}

    output_type operator()() const
    {
        return sum_;
    }

    HALMD_GPU_ENABLED void operator()(typename iterator::value_type const& tuple)
    {
        input_type value = halmd::get<0>(tuple);
        input_type factor = halmd::get<1>(tuple);
        sum_ += value * factor;
    }

    HALMD_GPU_ENABLED void operator()(sum_with_constant const& acc)
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
    typedef halmd::zip_iterator<input_type const*, input_type const*> iterator;

    // non-DefaultConstructible
    explicit sum_of_squares(output_type const& init) : sum_(init) {}

    output_type operator()() const
    {
        return sum_;
    }

    HALMD_GPU_ENABLED void operator()(typename iterator::value_type const& value)
    {
        input_type first = halmd::get<0>(value);
        input_type second = halmd::get<1>(value);
        sum_ += output_type(first) * output_type(second);
    }

    HALMD_GPU_ENABLED void operator()(sum_of_squares const& acc)
    {
        sum_ += acc.sum_;
    }

private:
    output_type sum_;
};

template <typename input_type, typename output_type>
class sum_of_squares_with_constant
{
public:
    typedef halmd::zip_iterator<input_type const*, input_type const*, halmd::constant_iterator<input_type> > iterator;

    // non-DefaultConstructible
    explicit sum_of_squares_with_constant(output_type const& init) : sum_(init) {}

    output_type operator()() const
    {
        return sum_;
    }

    HALMD_GPU_ENABLED void operator()(typename iterator::value_type const& value)
    {
        input_type first = halmd::get<0>(value);
        input_type second = halmd::get<1>(value);
        input_type factor = halmd::get<2>(value);
        sum_ += output_type(first) * output_type(second * factor);
    }

    HALMD_GPU_ENABLED void operator()(sum_of_squares_with_constant const& acc)
    {
        sum_ += acc.sum_;
    }

private:
    output_type sum_;
};

template <typename input_type, typename output_type>
class sum_of_cubes
{
public:
    typedef halmd::zip_iterator<input_type const*, input_type const*, input_type const*> iterator;

    // non-DefaultConstructible
    explicit sum_of_cubes(output_type const& init) : sum_(init) {}

    output_type operator()() const
    {
        return sum_;
    }

    HALMD_GPU_ENABLED void operator()(typename iterator::value_type const& value)
    {
        input_type first = halmd::get<0>(value);
        input_type second = halmd::get<1>(value);
        input_type third = halmd::get<2>(value);
        sum_ += output_type(first) * output_type(second) * output_type(third);
    }

    HALMD_GPU_ENABLED void operator()(sum_of_cubes const& acc)
    {
        sum_ += acc.sum_;
    }

private:
    output_type sum_;
};

#endif /* ! TEST_UNIT_ALGORITHM_GPU_REDUCE_KERNEL_HPP */
