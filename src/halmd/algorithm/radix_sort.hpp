/* Parallel radix sort
 *
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

#ifndef HALMD_ALGORITHM_RADIX_SORT_HPP
#define HALMD_ALGORITHM_RADIX_SORT_HPP

#include <algorithm>
#include <boost/ref.hpp>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/algorithm/prefix_sum.hpp>

namespace halmd
{

/*
 * Parallel radix sort
 */
template <typename T>
class radix_sort
{
public:
    typedef cuda::vector<unsigned int> key_vector;
    typedef cuda::vector<T> val_vector;

    enum { BUCKETS_PER_THREAD = gpu::radix_sort::BUCKETS_PER_THREAD };
    enum { BUCKET_SIZE = gpu::radix_sort::BUCKET_SIZE };
    enum { RADIX = gpu::radix_sort::RADIX };

public:
    radix_sort() : count_(0) {}

    /**
     * reallocate parallel radix sort for given element count
     */
    void resize(unsigned int count, unsigned int threads)
    {
        // don't reallocate the same count
        if (count_ == count) return;
        count_ = count;

        // compute optimal CUDA thread configuration
        int dev = cuda::device::get();
        unsigned int max_threads, max_blocks;
        cuda::device::properties prop(dev);
        max_threads = prop.multi_processor_count() * prop.max_threads_per_block();
        max_blocks = max_threads / (threads * BUCKETS_PER_THREAD / 2);
        blocks_ = std::min((count + 2 * threads - 1) / (2 * threads), max_blocks);
        threads_ = threads;

        // allocate global device memory
        g_scan.resize(blocks_ * threads_ * BUCKETS_PER_THREAD, BUCKET_SIZE);
        g_bucket.resize(blocks_ * threads_ * BUCKETS_PER_THREAD);
        g_key.resize(count);
        g_val.resize(count);
    }

    /**
     * radix sort given keys and values in-place
     */
    void operator()(key_vector& g_key_, val_vector& g_val_)
    {
        typedef boost::reference_wrapper<key_vector> key_ref;
        typedef boost::reference_wrapper<val_vector> val_ref;
        size_t const shmem = threads_ * BUCKETS_PER_THREAD * sizeof(unsigned int);

        assert(g_key_.size() == g_key.size());
        assert(g_val_.size() == g_val.size());

        // assign GPU dual buffers, as in the CUDA SDK radix sort example
        key_ref key[2] = { boost::ref(g_key_), boost::ref(g_key) };
        val_ref val[2] = { boost::ref(g_val_), boost::ref(g_val) };

        for (unsigned int r = 0; r < 32; r += RADIX) {
            key_vector const& ikey = key[0];
            key_vector& okey = key[1];
            val_vector const& ival = val[0];
            val_vector& oval = val[1];

            // compute partial radix counts
            cuda::configure(blocks_, threads_, shmem);
            gpu::radix_sort::histogram_keys(ikey, g_bucket, g_key.size(), r);
            // parallel prefix sum over radix counts
            g_scan(g_bucket);
            // permute array
            cuda::configure(blocks_, threads_, shmem);
            gpu::radix_sort::permute(ikey, okey, ival, oval, g_bucket, g_key.size(), r);

            // swap GPU dual buffers
            std::swap(key[0], key[1]);
            std::swap(val[0], val[1]);
        }
    }

private:
    unsigned int count_, blocks_, threads_;
    prefix_sum<unsigned int> g_scan;
    key_vector g_bucket, g_key;
    val_vector g_val;
};

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_RADIX_SORT_HPP */
