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
#ifndef HALMD_UTILITY_GPU_DSFLOAT_AS_FLOAT_HPP
#define HALMD_UTILITY_GPU_DSFLOAT_AS_FLOAT_HPP

#include <halmd/numeric/mp/dsfloat.hpp>

namespace halmd {

/*
 * Temporary wrapper class for obtaining a float4 pointer from cuda::vector and dsfloat_vector
 * in a unified manner. The dsfloat_vector version only uses the significant half of the dsfloats.
 * At a later point (after all kernels have been adjusted for dsfloat data types) an implicit conversion
 * between dsfloat_ptr and float4* could be recreated to avoid the need for this wrapper.
 */

namespace detail {

template<typename... Ts> struct make_void { typedef void type;};
template<typename... Ts> using void_t = typename make_void<Ts...>::type;

template<typename T, typename = void>
struct dsfloat_as_float_helper {
    static auto remove(T& vec) -> decltype(vec) {
        return  vec;
    }
};

template<typename T>
struct dsfloat_as_float_helper<T, void_t<decltype(std::declval<T const&>().storage())>> {
    static auto remove(T& vec) -> decltype(vec.storage()) {
        return  vec.storage();
    }
};

} // namespace detail

template<typename T>
auto dsfloat_as_float(T &vec) -> decltype(detail::dsfloat_as_float_helper<T>::remove(vec)) {
    return detail::dsfloat_as_float_helper<T>::remove(vec);
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_GPU_DSFLOAT_AS_FLOAT_HPP */