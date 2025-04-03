/*
 * Copyright © 2025  Felix Höfling
 * Copyright © 2017  Daniel Kirchner
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

#ifndef HALMD_NUMERIC_MP_DSFLOAT_VECTOR_HPP
#define HALMD_NUMERIC_MP_DSFLOAT_VECTOR_HPP

#include <halmd/numeric/mp/dsfloat.hpp>

namespace halmd {

/*
 * The template parameter must model a ContiguousContainer, e.g., std::vector.
 */
template<typename container_type>
class dsfloat_vector
{
public:
    typedef dsfloat_vector<container_type> vector_type;
    typedef typename container_type::value_type value_type;
    typedef dsfloat_ptr<value_type> pointer;
    typedef dsfloat_const_ptr<value_type> const const_pointer;
    typedef size_t size_type;

    dsfloat_vector(size_type size) : data_(size)
    {
        data_.reserve(size * 2);
    }

    size_type size() const
    {
        return data_.size();
    }

    void resize(size_type size)
    {
        data_.reserve(size * 2);
        data_.resize(size);
    }

    void swap(vector_type& v)
    {
        using std::swap;
        swap(data_, v.data_);
    }

    pointer data()
    {
        return pointer {
            &*data_.begin()
          , &*(data_.begin() + data_.size())
        };
    }

    operator pointer()
    {
        return data();
    }

    const_pointer data() const
    {
        return const_pointer {
            &*data_.begin()
          , &*(data_.begin() + data_.size())
        };
    }

    operator const_pointer() const
    {
        return data();
    }

    operator container_type const&() const
    {
        return data_;
    }

    operator container_type&()
    {
        return data_;
    }

private:
    container_type data_;
};

template<typename T>
inline void swap(dsfloat_vector<T>& lhs, dsfloat_vector<T>& rhs) noexcept
{
    lhs.swap(rhs);
}

} // namespace halmd

#endif /* ! HALMD_NUMERIC_MP_DSFLOAT_VECTOR_HPP */
