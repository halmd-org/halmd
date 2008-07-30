/* Parallel radix sort
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

#ifndef MDSIM_RADIX_HPP
#define MDSIM_RADIX_HPP

#include <algorithm>
#include <boost/ref.hpp>
#include "gpu/radix_glue.hpp"
#include "scan.hpp"

namespace mdsim
{

/*
 * Parallel radix sort
 */
template <typename T>
class radix_sort
{
public:
    /**
     * allocate parallel radix sort for given element count
     */
    radix_sort(uint const& count, uint const& blocks, uint const& threads)
    : count(count), blocks(blocks), threads(threads), scan(blocks * threads * gpu::radix::BUCKETS_PER_THREAD, gpu::radix::BUCKET_SIZE)
    {
	// allocate radix count buckets
	g_bucket.resize(blocks * threads * gpu::radix::BUCKETS_PER_THREAD);

	// allocate temporary sort key and value buffers
	g_key.resize(count);
	g_val.resize(count);
    }

    /**
     * radix sort given keys and values in-place
     */
    void operator()(cuda::vector<uint>& g_key_, cuda::vector<T>& g_val_, cuda::stream& stream)
    {
	assert(g_key_.size() == count);
	assert(g_val_.size() == count);

	using namespace boost;

	// assign GPU dual buffers, as in the CUDA SDK radix sort example
	reference_wrapper<cuda::vector<uint> > key[2] = { ref(g_key_), ref(g_key) };
	reference_wrapper<cuda::vector<T> > val[2] = { ref(g_val_), ref(g_val) };

	for (uint r = 0; r < 32; r += gpu::radix::RADIX) {
	    // compute partial radix counts
	    cuda::configure(blocks, threads, threads * gpu::radix::BUCKETS_PER_THREAD * sizeof(uint), stream);
	    gpu::radix::histogram_keys(key[0].get(), g_bucket, count, r);

	    // parallel prefix sum over radix counts
	    scan(g_bucket, stream);

	    // permute array
	    cuda::configure(blocks, threads, threads * gpu::radix::BUCKETS_PER_THREAD * sizeof(uint), stream);
	    gpu::radix::permute(key[0].get(), key[1].get(), val[0].get(), val[1].get(), g_bucket, count, r);

	    // swap GPU dual buffers
	    std::swap(key[0], key[1]);
	    std::swap(val[0], val[1]);
	}
    }

private:
    const uint count, blocks, threads;
    cuda::vector<uint> g_bucket, g_key;
    cuda::vector<T> g_val;
    prefix_sum<uint> scan;
};

} // namespace mdsim

#endif /* ! MDSIM_RADIX_HPP */
