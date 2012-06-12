/*
 * Copyright © 2012 Peter Colberg
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

#ifndef HALMD_UTILITY_ITERATOR_HPP
#define HALMD_UTILITY_ITERATOR_HPP

#include <halmd/config.hpp>
#ifdef __CUDACC__
# include <halmd/utility/tuple.hpp>
#endif

#include <iterator>
#ifndef HALMD_NO_CXX11
# include <tuple>
# include <type_traits>
#endif

namespace halmd {

/**
 * Constant value iterator.
 *
 * This iterator implements a partial Random Access Iterator that returns
 * a constant value when accessed at any given offset. In particular, the
 * iterator does not hold a distance, which avoids a duplicate distance
 * value when combined with a Random Access Iterator using a zip iterator.
 */
template <typename T>
class constant_iterator
{
public:
    typedef T const value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef typename std::iterator_traits<pointer>::difference_type difference_type;
    typedef typename std::iterator_traits<pointer>::iterator_category iterator_category;

    /**
     * Initialise with given value, or default construct.
     */
    constant_iterator(value_type const& value = value_type()) : value_(value) {}

#ifdef __CUDACC__
    /**
     * Returns reference to constant value.
     *
     * This function is invoked on the GPU.
     */
    HALMD_GPU_ENABLED reference operator[](difference_type const&) const
    {
        return value_;
    }
#endif

private:
    value_type value_;
};

/**
 * Zip iterator.
 *
 * This iterator combines multiple Random Access Iterators into a partial
 * Random Access Iterator. The iterator bridges the gap between the host,
 * where a std::tuple of iterators is used for construction, and the GPU,
 * where a halmd::tuple of value references is returned upon dereference
 * at an offset.
 */
template <typename first_iterator, typename second_iterator, typename third_iterator = void>
class zip_iterator;

template <typename first_iterator, typename second_iterator>
class zip_iterator<first_iterator, second_iterator>
{
private:
    typedef typename std::iterator_traits<first_iterator>::reference first_reference;
    typedef typename std::iterator_traits<second_iterator>::reference second_reference;

    struct enabler {}; // a private type avoids misuse

public:
    typedef typename std::iterator_traits<first_iterator>::difference_type difference_type;
    typedef typename std::iterator_traits<first_iterator>::iterator_category iterator_category;
    typedef tuple<first_reference, second_reference> value_type;
    typedef value_type* pointer;
    typedef value_type& reference;

    /**
     * Default constructor.
     */
    zip_iterator() {}

#ifndef HALMD_NO_CXX11
    /**
     * Initialise contained iterators from given tuple of iterators.
     */
    template <typename Iterator1, typename Iterator2>
    zip_iterator(std::tuple<Iterator1, Iterator2> const& iterator
      , typename std::enable_if<
            std::is_convertible<Iterator1, first_iterator>::value
            && std::is_convertible<Iterator2, second_iterator>::value
          , enabler
        >::type = enabler()
    )
      : first_(std::get<0>(iterator))
      , second_(std::get<1>(iterator))
    {}

    /**
     * Initialise first iterator from given iterator.
     *
     * This function is useful for construction of a last iterator,
     * which requires only the first iterator for distance calculation.
     */
    template <typename Iterator>
    zip_iterator(std::tuple<Iterator> const& iterator
      , typename std::enable_if<
            std::is_convertible<Iterator, first_iterator>::value
          , enabler
        >::type = enabler()
    )
      : first_(std::get<0>(iterator))
    {}
#endif

    /**
     * Returns distance to given iterator.
     *
     * This function is invoked on the host, and the result may be passed
     * to a GPU kernel function. This is maximally efficient, since it
     * avoids an unneeded copy of the second and subsequent iterators.
     */
    difference_type operator-(zip_iterator const& other) const
    {
        return first_ - other.first_;
    }

#ifdef __CUDACC__
    /**
     * Returns tuple of value references.
     *
     * This function is invoked on the GPU.
     */
    HALMD_GPU_ENABLED value_type operator[](difference_type const& offset) const
    {
        return tie(first_[offset], second_[offset]);
    }
#endif

private:
    first_iterator first_;
    second_iterator second_;
};

template <typename first_iterator, typename second_iterator, typename third_iterator>
class zip_iterator
{
private:
    typedef typename std::iterator_traits<first_iterator>::reference first_reference;
    typedef typename std::iterator_traits<second_iterator>::reference second_reference;
    typedef typename std::iterator_traits<third_iterator>::reference third_reference;

    struct enabler {}; // a private type avoids misuse

public:
    typedef typename std::iterator_traits<first_iterator>::difference_type difference_type;
    typedef typename std::iterator_traits<first_iterator>::iterator_category iterator_category;
    typedef tuple<first_reference, second_reference, third_reference> value_type;
    typedef value_type* pointer;
    typedef value_type& reference;

    /**
     * Default constructor.
     */
    zip_iterator() {}

#ifndef HALMD_NO_CXX11
    /**
     * Initialise contained iterators from given tuple of iterators.
     */
    template <typename Iterator1, typename Iterator2, typename Iterator3>
    zip_iterator(std::tuple<Iterator1, Iterator2, Iterator3> const& iterator
      , typename std::enable_if<
            std::is_convertible<Iterator1, first_iterator>::value
            && std::is_convertible<Iterator2, second_iterator>::value
            && std::is_convertible<Iterator3, third_iterator>::value
          , enabler
        >::type = enabler()
    )
      : first_(std::get<0>(iterator))
      , second_(std::get<1>(iterator))
      , third_(std::get<2>(iterator))
    {}

    /**
     * Initialise first iterator from given iterator.
     *
     * This function is useful for construction of a last iterator,
     * which requires only the first iterator for distance calculation.
     */
    template <typename Iterator>
    zip_iterator(std::tuple<Iterator> const& iterator
      , typename std::enable_if<
            std::is_convertible<Iterator, first_iterator>::value
          , enabler
        >::type = enabler()
    )
      : first_(std::get<0>(iterator))
    {}
#endif

    /**
     * Returns distance to given iterator.
     *
     * This function is invoked on the host, and the result may be passed
     * to a GPU kernel function. This is maximally efficient, since it
     * avoids an unneeded copy of the second and subsequent iterators.
     */
    difference_type operator-(zip_iterator const& other) const
    {
        return first_ - other.first_;
    }

#ifdef __CUDACC__
    /**
     * Returns tuple of value references.
     *
     * This function is invoked on the GPU.
     */
    HALMD_GPU_ENABLED value_type operator[](difference_type const& offset) const
    {
        return tie(first_[offset], second_[offset], third_[offset]);
    }
#endif

private:
    first_iterator first_;
    second_iterator second_;
    third_iterator third_;
};

} // namespace halmd

#endif /* ! HALMD_UTILITY_ITERATOR_HPP */
