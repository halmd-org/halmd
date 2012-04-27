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

#ifndef TEST_UNIT_MDSIM_POSITIONS_LATTICE_ITERATOR_HPP
#define TEST_UNIT_MDSIM_POSITIONS_LATTICE_ITERATOR_HPP

#include <cstddef> // std::size_t, std::ptrdiff_t
#include <iterator>

#include <halmd/mdsim/positions/lattice_primitive.hpp>

/**
 * Generate square/cubic lattice vectors with unit lattice constant.
 */
template <typename vector_type, typename primitive_type = halmd::sc_lattice_primitive>
class lattice_iterator
{
public:
    typedef std::forward_iterator_tag iterator_category;
    typedef std::ptrdiff_t difference_type;
    typedef vector_type* pointer;
    typedef vector_type& reference;
    typedef vector_type value_type;

    /**
     * Increment lattice site index.
     */
    lattice_iterator& operator++()
    {
        ++index_;
        return *this;
    }

    /**
     * Generate lattice vector.
     *
     * The components are multiples of integers, with each component
     * shifted by +0.5 to the centre of the unit cell. Therefore the
     * lattice vector is exactly representable in floating-point.
     */
    vector_type operator*()
    {
        vector_type r;
        primitive_type()(r, nsite_, index_);
        return r;
    }

    bool operator!=(lattice_iterator const& other)
    {
        return index_ != other.index_;
    }

    lattice_iterator(std::size_t nparticle, std::size_t index) : nsite_(nsite(nparticle)), index_(index) {}

private:
    /**
     * Convert number of particles to number of lattice sites per dimension.
     */
    static std::size_t nsite(std::size_t nparticle)
    {
        return ceil(pow(nparticle, 1. / vector_type::static_size));
    }

    /** number of lattice sites per dimension */
    std::size_t nsite_;
    /** linear index of lattice site */
    std::size_t index_;
};

#endif /* TEST_UNIT_MDSIM_POSITIONS_LATTICE_ITERATOR_HPP */
