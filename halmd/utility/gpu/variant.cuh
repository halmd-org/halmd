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

#ifndef HALMD_UTILITY_GPU_VARIANT_CUH
#define HALMD_UTILITY_GPU_VARIANT_CUH

#include <boost/mpl/at.hpp>
#include <boost/mpl/begin.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/next.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/value_type.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

namespace halmd
{
namespace utility { namespace gpu
{

// import into current namespace
using boost::mpl::set;
using boost::mpl::map;
using boost::mpl::pair;
using boost::mpl::int_;

namespace detail
{

template <
    typename Sequence
  , typename Begin
  , typename End
  , typename Enable = void
>
union variant_iterate_range {};

template <
    typename Sequence
  , typename Begin
  , typename End
>
union variant_iterate_range<
    Sequence
  , Begin
  , End
  , typename boost::disable_if<boost::is_same<Begin, End> >::type
>
{
    typename boost::mpl::value_type<
        Sequence
      , typename boost::mpl::deref<Begin>::type
    >::type head;
    variant_iterate_range<
        Sequence
      , typename boost::mpl::next<Begin>::type
      , End
    > range;
};

} // namespace detail

template <typename Sequence>
union variant
{
    detail::variant_iterate_range<
        Sequence
      , typename boost::mpl::begin<Sequence>::type
      , typename boost::mpl::end<Sequence>::type
    > range;
};

template <typename Key, typename Sequence>
__device__ __host__ typename boost::enable_if<
    boost::mpl::has_key<Sequence, Key>
  , typename boost::mpl::at<Sequence, Key>::type const&
>::type
inline get(variant<Sequence> const& value)
{
    return reinterpret_cast<
        typename boost::mpl::at<Sequence, Key>::type const&
    >(value);
}

template <int Key, typename Sequence>
__device__ __host__ typename boost::enable_if<
    boost::mpl::has_key<Sequence, boost::mpl::int_<Key> >
  , typename boost::mpl::at<Sequence, boost::mpl::int_<Key> >::type const&
>::type
inline get(variant<Sequence> const& value)
{
    return reinterpret_cast<
        typename boost::mpl::at<Sequence, boost::mpl::int_<Key> >::type const&
    >(value);
}

template <typename Key, typename Sequence, int dim, enum cudaTextureReadMode mode>
__device__ __host__ typename boost::enable_if<
    boost::mpl::has_key<Sequence, Key>
  , texture<
        typename boost::mpl::at<Sequence, Key>::type
      , dim
      , mode
    > const&
>::type
inline get(texture<variant<Sequence>, dim, mode> const& value)
{
    return reinterpret_cast<
        texture<
            typename boost::mpl::at<Sequence, Key>::type
          , dim
          , mode
        > const&
    >(value);
}

template <int Key, typename Sequence, int dim, enum cudaTextureReadMode mode>
__device__ __host__ typename boost::enable_if<
    boost::mpl::has_key<Sequence, boost::mpl::int_<Key> >
  , texture<
        typename boost::mpl::at<Sequence, boost::mpl::int_<Key> >::type
      , dim
      , mode
    > const&
>::type
inline get(texture<variant<Sequence>, dim, mode> const& value)
{
    return reinterpret_cast<
        texture<
            typename boost::mpl::at<Sequence, boost::mpl::int_<Key> >::type
          , dim
          , mode
        > const&
    >(value);
}

}} // namespace utility::gpu

} // namespace halmd

#endif /* ! HALMD_UTILITY_GPU_VARIANT_CUH */
