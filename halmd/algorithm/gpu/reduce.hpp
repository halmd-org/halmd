/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_REDUCE_HPP
#define HALMD_ALGORITHM_GPU_REDUCE_HPP

#include <algorithm>
#include <numeric>
#include <stdexcept>

#include <boost/lambda/lambda.hpp> // before other Boost.Lambda headers
#include <boost/lambda/casts.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/algorithm/gpu/reduce_kernel.hpp>

namespace halmd
{
namespace algorithm { namespace gpu
{

template <
    typename reduce_transform
  , typename host_output_type
>
struct reduce_blocks;

template <
    typename reduce_transform
  , typename input_type
  , typename coalesced_input_type       = input_type
  , typename output_type                = input_type
  , typename coalesced_output_type      = output_type
  , typename host_output_type           = coalesced_output_type
  , typename input_transform            = identity_
  , typename output_transform           = identity_
>
struct reduce
{
    typedef reduce_wrapper<
        reduce_transform
      , input_type
      , coalesced_input_type
      , output_type
      , coalesced_output_type
      , input_transform
      , output_transform
    > wrapper_type;

    typedef cuda::function<
        void (coalesced_input_type const*, coalesced_output_type*, unsigned int)
    > reduce_impl_type;

    reduce(int blocks = 16, int threads = (64 << DEVICE_SCALE))
      : dim(blocks, threads)
      , g_block(blocks)
      , h_block(blocks)
      , reduce_impl(get_reduce_impl(threads))
    {}

    host_output_type operator()(cuda::vector<coalesced_input_type> const& g_v)
    {
        cuda::configure(dim.grid, dim.block);
        reduce_impl(g_v, g_block, g_v.size());
        cuda::copy(g_block, h_block);
        return reduce_blocks<reduce_transform, host_output_type>()(h_block);
    }

    static reduce_impl_type get_reduce_impl(int threads)
    {
        switch (threads) {
          case 512:
            return wrapper_type::kernel.reduce_impl_512;
          case 256:
            return wrapper_type::kernel.reduce_impl_256;
          case 128:
            return wrapper_type::kernel.reduce_impl_128;
          case 64:
            return wrapper_type::kernel.reduce_impl_64;
          case 32:
            return wrapper_type::kernel.reduce_impl_32;
          default:
            throw std::logic_error("invalid reduction thread count");
        }
    }

    cuda::config dim;
    cuda::vector<coalesced_output_type> g_block;
    cuda::host::vector<coalesced_output_type> h_block;
    reduce_impl_type reduce_impl;
};

template <typename host_output_type>
struct reduce_blocks<sum_, host_output_type>
{
    template <typename coalesced_output_type>
    host_output_type operator()(cuda::host::vector<coalesced_output_type> const& h_block)
    {
        using namespace boost::lambda;
        using boost::lambda::_1;
        return std::accumulate(
            boost::make_transform_iterator(
                h_block.begin()
              , ret<host_output_type>(ll_static_cast<host_output_type>(_1))
            )
          , boost::make_transform_iterator(
                h_block.end()
              , ret<host_output_type>(ll_static_cast<host_output_type>(_1))
            )
          , host_output_type(0)
        );
    }
};

template <typename host_output_type>
struct reduce_blocks<max_, host_output_type>
{
    template <typename coalesced_output_type>
    host_output_type operator()(cuda::host::vector<coalesced_output_type> const& h_block)
    {
        using namespace boost::lambda;
        using boost::lambda::_1;
        return *std::max_element(
            boost::make_transform_iterator(
                h_block.begin()
              , ret<host_output_type>(ll_static_cast<host_output_type>(_1))
            )
          , boost::make_transform_iterator(
                h_block.end()
              , ret<host_output_type>(ll_static_cast<host_output_type>(_1))
            )
        );
    }
};

}} // namespace algorithm::gpu

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_REDUCE_HPP */
