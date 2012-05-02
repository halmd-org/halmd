/*
 * Copyright © 2008, 2012  Peter Colberg
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

#ifndef HALMD_MDSIM_POSITIONS_LATTICE_PRIMITIVE_HPP
#define HALMD_MDSIM_POSITIONS_LATTICE_PRIMITIVE_HPP

#include <boost/utility/enable_if.hpp>

#include <halmd/config.hpp>
#include <halmd/utility/multi_index.hpp>

namespace halmd {

template <typename Position, typename Shape, typename Enable = void>
class close_packed_lattice;

/**
 * Two-dimensional hexagonal close-packed unit lattice.
 *
 * @tparam Position 2-dimensional floating-n array type
 * @tparam Shape 2-dimensional integer array type
 */
template <typename Position, typename Shape>
class close_packed_lattice<Position, Shape
  , typename boost::enable_if_c<(Position::static_size == 2)>::type>
{
public:
    typedef Position result_type;
    typedef Shape shape_type;
    typedef typename shape_type::size_type size_type;

    /**
     * Construct unit lattice primitive of given shape.
     *
     * @param shape number of unit cells per dimension
     */
    explicit HALMD_GPU_ENABLED close_packed_lattice(shape_type const& shape) : shape_(shape) {}

    /**
     * Returns lattice position for given particle index.
     *
     * @param n particle index with 0 ≤ index < size()
     */
    HALMD_GPU_ENABLED result_type operator()(size_type n) const
    {
        shape_type cell = offset_to_multi_index(n >> 1, shape_);
        result_type r;
        r[0] = cell[0] + ((n & 1) * 2 + 1) / float_type(4);
        r[1] = cell[1] + ((n & 1) * 2 + 1) / float_type(4);
        return r;
    }

    /**
     * Returns number of unit cells per dimension.
     */
    HALMD_GPU_ENABLED shape_type const& shape() const
    {
        return shape_;
    }

    /**
     * Returns total number of lattice points.
     */
    HALMD_GPU_ENABLED size_type size() const
    {
        return 2 * shape_[0] * shape_[1];
    }

private:
    typedef typename result_type::value_type float_type;

    /** number of unit cells per dimension */
    shape_type shape_;
};

/**
 * Three-dimensional face centered cubic unit lattice.
 *
 * @tparam Position 3-dimensional floating-n array type
 * @tparam Shape 3-dimensional integer array type
 */
template <typename Position, typename Shape>
class close_packed_lattice<Position, Shape
  , typename boost::enable_if_c<(Position::static_size == 3)>::type>
{
public:
    typedef Position result_type;
    typedef Shape shape_type;
    typedef typename shape_type::size_type size_type;

    /**
     * Construct unit lattice primitive of given shape.
     *
     * @param shape number of unit cells per dimension
     */
    explicit HALMD_GPU_ENABLED close_packed_lattice(shape_type const& shape) : shape_(shape) {}

    /**
     * Returns lattice position for given particle index.
     *
     * @param n particle index with 0 ≤ index < size()
     */
    HALMD_GPU_ENABLED result_type operator()(size_type n) const
    {
        shape_type cell = offset_to_multi_index(n >> 2, shape_);
        result_type r;
        r[0] = cell[0] + (((n ^ (n >> 1)) & 1) * 2 + 1) / float_type(4);
        r[1] = cell[1] + ((n & 1) * 2 + 1) / float_type(4);
        r[2] = cell[2] + ((n & 2) + 1) / float_type(4);
        return r;
    }

    /**
     * Returns number of unit cells per dimension.
     */
    HALMD_GPU_ENABLED shape_type const& shape() const
    {
        return shape_;
    }

    /**
     * Returns total number of lattice points.
     */
    HALMD_GPU_ENABLED size_type size() const
    {
        return 4 * shape_[0] * shape_[1] * shape_[2];
    }

private:
    typedef typename result_type::value_type float_type;

    /** number of unit cells per dimension */
    shape_type shape_;
};

template <typename Position, typename Shape, typename Enable = void>
class primitive_lattice;

/**
 * Two-dimensional square unit lattice.
 *
 * @tparam Position 2-dimensional floating-n array type
 * @tparam Shape 2-dimensional integer array type
 */
template <typename Position, typename Shape>
class primitive_lattice<Position, Shape
  , typename boost::enable_if_c<(Position::static_size == 2)>::type>
{
public:
    typedef Position result_type;
    typedef Shape shape_type;
    typedef typename shape_type::size_type size_type;

    /**
     * Construct unit lattice primitive of given shape.
     *
     * @param shape number of unit cells per dimension
     */
    explicit HALMD_GPU_ENABLED primitive_lattice(shape_type const& shape) : shape_(shape) {}

    /**
     * Returns lattice position for given particle index.
     *
     * @param n particle index with 0 ≤ index < size()
     */
    HALMD_GPU_ENABLED result_type operator()(size_type n) const
    {
        shape_type cell = offset_to_multi_index(n, shape_);
        result_type r;
        r[0] = cell[0] + float_type(0.5);
        r[1] = cell[1] + float_type(0.5);
        return r;
    }

    /**
     * Returns number of unit cells per dimension.
     */
    HALMD_GPU_ENABLED shape_type const& shape() const
    {
        return shape_;
    }

    /**
     * Returns total number of lattice points.
     */
    HALMD_GPU_ENABLED size_type size() const
    {
        return shape_[0] * shape_[1];
    }

private:
    typedef typename result_type::value_type float_type;

    /** number of unit cells per dimension */
    shape_type shape_;
};

/**
 * Three-dimensional primitive cubic unit lattice.
 *
 * @tparam Position 3-dimensional floating-n array type
 * @tparam Shape 3-dimensional integer array type
 */
template <typename Position, typename Shape>
class primitive_lattice<Position, Shape
  , typename boost::enable_if_c<(Position::static_size == 3)>::type>
{
public:
    typedef Position result_type;
    typedef Shape shape_type;
    typedef typename shape_type::size_type size_type;

    /**
     * Construct unit lattice primitive of given shape.
     *
     * @param shape number of unit cells per dimension
     */
    explicit HALMD_GPU_ENABLED primitive_lattice(shape_type const& shape) : shape_(shape) {}

    /**
     * Returns lattice position for given particle index.
     *
     * @param n particle index with 0 ≤ index < size()
     */
    HALMD_GPU_ENABLED result_type operator()(size_type n) const
    {
        shape_type cell = offset_to_multi_index(n, shape_);
        result_type r;
        r[0] = cell[0] + float_type(0.5);
        r[1] = cell[1] + float_type(0.5);
        r[2] = cell[2] + float_type(0.5);
        return r;
    }

    /**
     * Returns number of unit cells per dimension.
     */
    HALMD_GPU_ENABLED shape_type const& shape() const
    {
        return shape_;
    }

    /**
     * Returns total number of lattice points.
     */
    HALMD_GPU_ENABLED size_type size() const
    {
        return shape_[0] * shape_[1] * shape_[2];
    }

private:
    typedef typename result_type::value_type float_type;

    /** number of unit cells per dimension */
    shape_type shape_;
};

/**
 * Four-dimensional primitive tesseractic unit lattice.
 *
 * @tparam Position 4-dimensional floating-n array type
 * @tparam Shape 4-dimensional integer array type
 */
template <typename Position, typename Shape>
class primitive_lattice<Position, Shape
  , typename boost::enable_if_c<(Position::static_size == 4)>::type>
{
public:
    typedef Position result_type;
    typedef Shape shape_type;
    typedef typename shape_type::size_type size_type;

    /**
     * Construct unit lattice primitive of given shape.
     *
     * @param shape number of unit cells per dimension
     */
    explicit HALMD_GPU_ENABLED primitive_lattice(shape_type const& shape) : shape_(shape) {}

    /**
     * Returns lattice position for given particle index.
     *
     * @param n particle index with 0 ≤ index < size()
     */
    HALMD_GPU_ENABLED result_type operator()(size_type n) const
    {
        shape_type cell = offset_to_multi_index(n, shape_);
        result_type r;
        r[0] = cell[0] + float_type(0.5);
        r[1] = cell[1] + float_type(0.5);
        r[2] = cell[2] + float_type(0.5);
        r[3] = cell[3] + float_type(0.5);
        return r;
    }

    /**
     * Returns number of unit cells per dimension.
     */
    HALMD_GPU_ENABLED shape_type const& shape() const
    {
        return shape_;
    }

    /**
     * Returns total number of lattice points.
     */
    HALMD_GPU_ENABLED size_type size() const
    {
        return shape_[0] * shape_[1] * shape_[2] * shape_[3];
    }

private:
    typedef typename result_type::value_type float_type;

    /** number of unit cells per dimension */
    shape_type shape_;
};

} // namespace halmd

#endif /* ! HALMD_MDSIM_POSITIONS_LATTICE_PRIMITIVE_HPP */
