/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_TUPLE_CUH
#define HALMD_ALGORITHM_GPU_TUPLE_CUH

#include <boost/mpl/int.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

//
// Poor man's tuple library for CUDA.
//
// This library tries to partially model Boost.Tuple using a very
// simple implementation. We cannot use Boost.Tuple because nvcc
// refuses to call perfectly valid (and most probably supported)
// functions if they do not have the __device__ attribute.
// Sigh…
//

namespace halmd
{
namespace algorithm { namespace gpu
{

// forward declaration
template <typename T0, typename T1 = void>
struct tuple;

template <int i, typename T0>
typename boost::enable_if<boost::is_same<boost::mpl::int_<0>, boost::mpl::int_<i> >, T0>::type
__device__ get(tuple<T0> const& t)
{
    return t.t0;
}

template <int i, typename T0, typename T1>
typename boost::enable_if<boost::is_same<boost::mpl::int_<0>, boost::mpl::int_<i> >, T0>::type
__device__ get(tuple<T0, T1> const& t)
{
    return t.t0;
}

template <int i, typename T0, typename T1>
typename boost::enable_if<boost::is_same<boost::mpl::int_<1>, boost::mpl::int_<i> >, T1>::type
__device__ get(tuple<T0, T1> const& t)
{
    return t.t1;
}

template <typename T0>
struct tuple<T0>
{
    T0 t0;
    friend T0 get<0>(tuple<T0> const&);
public:
    __device__ tuple(T0 t0) : t0(t0) {}
    __device__ tuple() {}
    template <typename TT0>
    __device__ tuple(tuple<TT0> const& t)
      : t0(t.t0) {}
    template <typename TT0>
    __device__ tuple& operator=(tuple<TT0> const& t) {
        t0 = get<0>(t); return *this;
    }
};

template <typename T0, typename T1>
struct tuple
{
    T0 t0; T1 t1;
    friend T0 get<0>(tuple<T0, T1> const&);
    friend T1 get<1>(tuple<T0, T1> const&);
public:
    __device__ tuple(T0 t0, T1 t1) : t0(t0), t1(t1) {}
    __device__ tuple(T0 t0) : t0(t0) {}
    __device__ tuple() {}
    template <typename TT0, typename TT1>
    __device__ tuple(tuple<TT0, TT1> const& t)
      : t0(t.t0), t1(t.t1) {}
    template <typename TT0, typename TT1>
    __device__ tuple& operator=(tuple<TT0, TT1> const& t) {
        t0 = get<0>(t); t1 = get<1>(t); return *this;
    }
};

template <typename T0>
__device__ tuple<T0> make_tuple(T0 t0)
{
    return tuple<T0>(t0);
}

template <typename T0, typename T1>
__device__ tuple<T0, T1> make_tuple(T0 t0, T1 t1)
{
    return tuple<T0, T1>(t0, t1);
}

template <typename T0>
__device__ tuple<T0&> tie(T0& t0)
{
    return tuple<T0&>(t0);
}

template <typename T0, typename T1>
__device__ tuple<T0&, T1&> tie(T0& t0, T1& t1)
{
    return tuple<T0&, T1&>(t0, t1);
}

}} // namespace algorithm::gpu

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_TUPLE_CUH */
