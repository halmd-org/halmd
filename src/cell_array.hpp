/* Molecular Dynamics simulation cell array
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

#ifndef MDSIM_CELL_ARRAY_HPP
#define MDSIM_CELL_ARRAY_HPP

#include <boost/multi_array.hpp>
#include "vector2d.hpp"
#include "vector3d.hpp"


namespace mdsim
{

template <typename T0, typename T1>
class cell_array;


/**
 * 2-dimensional Molecular Dynamics simulation cell array
 */
template <typename T0, typename T1>
class cell_array<T0, vector2d<T1> > : public boost::multi_array<T0, 2>
{
public:
    typedef T0* iterator;
    typedef T0 const* const_iterator;

protected:
    typedef boost::multi_array<T0, 2> _Base;

public:
    cell_array() : _Base() {}

    cell_array(size_t const& ncell) : _Base(boost::extents[ncell][ncell][ncell]) {}

    iterator begin()
    {
	return iterator(_Base::data());
    }

    const_iterator begin() const
    {
	return const_iterator(_Base::data());
    }

    iterator end()
    {
	return iterator(_Base::data() + _Base::num_elements());
    }

    const_iterator end() const
    {
	return const_iterator(_Base::data() + _Base::num_elements());
    }

    void resize(size_t ncell)
    {
	_Base::resize(boost::extents[ncell][ncell]);
    }

    T0& operator()(vector2d<T1> const& v)
    {
	return (*this)[size_t(v.x)][size_t(v.y)];
    }
};


/**
 * 3-dimensional Molecular Dynamics simulation cell array
 */
template <typename T0, typename T1>
class cell_array<T0, vector3d<T1> > : public boost::multi_array<T0, 3>
{
public:
    typedef T0* iterator;
    typedef T0 const* const_iterator;

protected:
    typedef boost::multi_array<T0, 3> _Base;

public:
    cell_array() : _Base() {}

    cell_array(size_t const& ncell) : _Base(boost::extents[ncell][ncell][ncell]) {}

    iterator begin()
    {
	return iterator(_Base::data());
    }

    const_iterator begin() const
    {
	return const_iterator(_Base::data());
    }

    iterator end()
    {
	return iterator(_Base::data() + _Base::num_elements());
    }

    const_iterator end() const
    {
	return const_iterator(_Base::data() + _Base::num_elements());
    }

    void resize(size_t ncell)
    {
	_Base::resize(boost::extents[ncell][ncell][ncell]);
    }

    T0& operator()(vector3d<T1> const& v)
    {
	return (*this)[size_t(v.x)][size_t(v.y)][size_t(v.z)];
    }
};

}

#endif /* ! MDSIM_CELL_ARRAY_HPP */
