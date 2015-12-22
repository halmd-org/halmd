/*
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_RADIX_SORT_HPP
#define HALMD_ALGORITHM_GPU_RADIX_SORT_HPP

#include <halmd/algorithm/gpu/radix_sort_kernel.hpp>
#include <halmd/algorithm/gpu/scan.hpp>

#include <algorithm>
#include <iterator>
#include <type_traits>

namespace halmd {
namespace algorithm {
namespace gpu {

/*
 * Parallel radix sort
 */
class radix_sort
{
public:
    /**
     * allocate parallel radix sort for given element count
     */
    radix_sort(unsigned int count, unsigned int threads)
      : count_(count)
      , threads_(threads)
      , blocks_(blocks(count_, threads_))
      , scan_(blocks_ * threads_ * BUCKETS_PER_THREAD, BUCKET_SIZE)
      , g_bucket_(blocks_ * threads_ * BUCKETS_PER_THREAD) {}

    /**
     * radix sort given keys in-place
     */
    template <typename Iterator>
    void operator()(
        Iterator const& first
      , Iterator const& last
    )
    {
        typedef typename std::iterator_traits<Iterator>::value_type key_type;
        std::size_t const shared_memory = threads_ * BUCKETS_PER_THREAD * sizeof(unsigned int);
        std::size_t const count = last - first;

        // do nothing in case of an empty array
        if (!count) return;

        // assign GPU dual buffers, as in the CUDA SDK radix sort example
        cuda::vector<key_type> g_key(count);
        std::pair<Iterator, Iterator> key = {first, g_key.begin()};

        for (unsigned int shift = 0; shift < 32; shift += RADIX) {
            // compute partial radix counts
            cuda::configure(blocks_, threads_, shared_memory);
            radix_sort_wrapper::kernel.histogram_key(
                &*key.first
              , g_bucket_
              , count
              , shift
            );
            // parallel prefix sum over radix counts
            scan_(g_bucket_);
            // permute array
            cuda::configure(blocks_, threads_, shared_memory);
            radix_sort_wrapper::kernel.permute_key(
                &*key.first
              , &*key.second
              , g_bucket_
              , count
              , shift
            );

            // swap GPU dual buffers
            std::swap(key.first, key.second);
        }
        cuda::thread::synchronize();
    }

    /**
     * radix sort given keys and values in-place
     */
    template <typename Iterator1, typename Iterator2>
    Iterator2 operator()(
        Iterator1 const& first1
      , Iterator1 const& last1
      , Iterator2 const& first2
    )
    {
        typedef typename std::iterator_traits<Iterator1>::value_type key_type;
        typedef typename std::iterator_traits<Iterator2>::value_type value_type;
        std::size_t const shared_memory = threads_ * BUCKETS_PER_THREAD * sizeof(unsigned int);
        std::size_t const count = last1 - first1;

        // assign GPU dual buffers, as in the CUDA SDK radix sort example
        cuda::vector<key_type> g_key(count);
        std::pair<Iterator1, Iterator1> key = {first1, g_key.begin()};
        cuda::vector<value_type> g_value(count);
        std::pair<Iterator2, Iterator2> value = {first2, g_value.begin()};

        // do nothing in case of an empty array
        if (!count) return value.first;

        for (unsigned int shift = 0; shift < 32; shift += RADIX) {
            // compute partial radix counts
            cuda::configure(blocks_, threads_, shared_memory);
            radix_sort_wrapper::kernel.histogram_key(
                &*key.first
              , g_bucket_
              , count
              , shift
            );
            // parallel prefix sum over radix counts
            scan_(g_bucket_);
            // permute array
            cuda::configure(blocks_, threads_, shared_memory);
            radix_sort_wrapper::kernel.permute_key_value(
                &*key.first
              , &*key.second
              , g_bucket_
              , count
              , shift
              , &*value.first
              , &*value.second
            );

            // swap GPU dual buffers
            std::swap(key.first, key.second);
            std::swap(value.first, value.second);
        }
        cuda::thread::synchronize();
        return value.first + count;
    }

private:
    /**
     * compute optimal CUDA execution configuration
     */
    static unsigned int blocks(unsigned int count, unsigned int threads)
    {
        int dev = cuda::device::get();
        unsigned int max_threads, max_blocks;
        cuda::device::properties prop(dev);
        max_threads = prop.multi_processor_count() * prop.max_threads_per_block();
        max_blocks = max_threads / (threads * BUCKETS_PER_THREAD / 2);
        return std::min((count + 2 * threads - 1) / (2 * threads), max_blocks);
    }

    unsigned int count_;
    unsigned int threads_;
    unsigned int blocks_;
    scan<unsigned int> scan_;
    cuda::vector<unsigned int> g_bucket_;
};

} // namespace gpu
} // namespace algorithm

/**
 * Radix sort keys in-place.
 */
template <typename Iterator>
inline typename std::enable_if<
    std::is_same<
        typename std::iterator_traits<Iterator>::value_type
      , unsigned int
    >::value
    && std::is_convertible<
        typename std::iterator_traits<Iterator>::iterator_category
      , cuda::device_random_access_iterator_tag
    >::value
  , void>::type
radix_sort(
    Iterator const& first
  , Iterator const& last
)
{
    unsigned int constexpr threads = 128;
    algorithm::gpu::radix_sort sort(last - first, threads);
    sort(first, last);
}

/**
 * Radix sort keys and values in-place.
 */
template <typename Iterator1, typename Iterator2>
inline typename std::enable_if<
    std::is_same<
        typename std::iterator_traits<Iterator1>::value_type
      , unsigned int
    >::value
    && std::is_same<
        typename std::iterator_traits<Iterator2>::value_type
      , unsigned int
    >::value
    && std::is_convertible<
        typename std::iterator_traits<Iterator1>::iterator_category
      , cuda::device_random_access_iterator_tag
    >::value
    && std::is_convertible<
        typename std::iterator_traits<Iterator2>::iterator_category
      , cuda::device_random_access_iterator_tag
    >::value
  , Iterator2>::type
radix_sort(
    Iterator1 const& first1
  , Iterator1 const& last1
  , Iterator2 const& first2
)
{
    unsigned int constexpr threads = 128;
    algorithm::gpu::radix_sort sort(last1 - first1, threads);
    return sort(first1, last1, first2);
}

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_RADIX_SORT_HPP */
