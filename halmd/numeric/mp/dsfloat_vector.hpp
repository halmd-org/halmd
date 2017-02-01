/*
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
namespace detail {
namespace numeric {
namespace mp {

template<typename T>
class dsfloat_vector;

} // namespace mp
} // namespace numeric
} // namespace detail

// import into top-level namespace
using detail::numeric::mp::dsfloat_vector;

namespace detail {
namespace numeric {
namespace mp {

template<typename T>
class dsfloat_vector {
public:
    typedef dsfloat_vector<T> vector_type;
    typedef T value_type;
    typedef dsfloat_ptr<T> pointer;
    typedef dsfloat_const_ptr<T> const const_pointer;
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
        data_.swap(v.data_);
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

    operator cuda::vector<T> const&() const {
        return data_;
    }

    operator cuda::vector<T>&() {
        return data_;
    }

private:
    cuda::vector<T> data_;
};

} // namespace mp
} // namespace numeric
} // namespace detail
} // namespace halmd

#endif /* ! HALMD_NUMERIC_MP_DSFLOAT_VECTOR_HPP */
