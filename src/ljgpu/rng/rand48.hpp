/* Parallelized rand48 random number generator for CUDA
 *
 * Copyright © 2007-2009  Peter Colberg
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

#ifndef LJGPU_RNG_RAND48_HPP
#define LJGPU_RNG_RAND48_HPP

#include <algorithm>
#include <iostream>
#include <ljgpu/algorithm/prefix_sum.hpp>
#include <ljgpu/rng/gpu/rand48.hpp>

namespace ljgpu
{

/**
 * parallelized rand48 random number generator for CUDA
 */
class rand48
{
public:
    /** type for saving or restoring generator state in memory */
    typedef ushort3 state_type;

public:
    rand48() {}

    /**
     * initialize random number generator with CUDA execution dimensions
     */
    rand48(cuda::config const& dim) : dim_(dim), state_(dim.threads()) {}

    /**
     * change random number generator CUDA execution dimensions
     */
    void resize(cuda::config const& dim)
    {
	if (state_.size() > 0) {
	    ushort3 x;
	    // save generator state using old dimensions
	    save(x);
	    // set new CUDA execution dimensions
	    dim_ = dim;
	    // reallocate global device memory for generator state
	    state_.resize(dim_.threads());
	    // restore generator state using new dimensions
	    restore(x);
	}
	else {
	    // set new CUDA execution dimensions
	    dim_ = dim;
	    // reallocate global device memory for generator state
	    state_.resize(dim_.threads());
	}
    }

    /**
     * initialize generator with 32-bit integer seed
     */
    void set(uint seed, cuda::stream& stream)
    {
	// compute leapfrog multipliers for initialization
	cuda::vector<uint48> g_a(dim_.threads()), g_c(dim_.threads());
	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::leapfrog(g_a);

	// compute leapfrog addends for initialization
	cuda::copy(g_a, g_c, stream);
	prefix_sum<uint48> scan(g_c.size(), dim_.threads_per_block());
	scan(g_c, stream);

	// initialize generator with seed
	cuda::vector<uint48> a(1), c(1);
	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::set(state_, g_a, g_c, a, c, seed);
	stream.synchronize();

	// copy leapfrog multiplier into constant device memory
	cuda::copy(a, gpu::rand48::a);
	// copy leapfrog addend into constant device memory
	cuda::copy(c, gpu::rand48::c);
    }

    /*
     * fill array with uniform random numbers in [0.0, 1.0)
     */
    void uniform(cuda::vector<float>& r, cuda::stream& stream)
    {
	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::uniform(state_, r, r.size());
    }

    /**
     * fill array with random integers in [0, 2^32-1]
     */
    void get(cuda::vector<uint>& r, cuda::stream& stream)
    {
	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::get(state_, r, r.size());
    }

    /**
     * generate random 2-dimensional Maxwell-Boltzmann distributed velocities
     */
    template <typename T>
    void boltzmann(cuda::vector<T>& v, float const& temperature, cuda::stream& stream)
    {
	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::boltzmann(v, temperature, state_);
    }

    /**
     * save generator state to memory
     */
    void save(state_type& mem)
    {
	cuda::stream stream;
	cuda::vector<ushort3> buf_gpu(1);
	cuda::host::vector<ushort3> buf(1);

	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::save(state_, buf_gpu);
	cuda::copy(buf_gpu, buf, stream);
	stream.synchronize();

	mem = buf[0];
    }

    /**
     * restore generator state from memory
     */
    void restore(state_type const& mem)
    {
	cuda::stream stream;

	// compute leapfrog multipliers for initialization
	cuda::vector<uint48> g_a(dim_.threads()), g_c(dim_.threads());
	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::leapfrog(g_a);

	// compute leapfrog addends for initialization
	cuda::copy(g_a, g_c, stream);
	prefix_sum<uint48> scan(g_c.size(), dim_.threads_per_block());
	scan(g_c, stream);

	// initialize generator from state
	cuda::vector<uint48> a(1), c(1);
	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::restore(state_, g_a, g_c, a, c, mem);
	stream.synchronize();

	// copy leapfrog multiplier into constant device memory
	cuda::copy(a, gpu::rand48::a);
	// copy leapfrog addend into constant device memory
	cuda::copy(c, gpu::rand48::c);
    }

    /**
     * save generator state to text-mode output stream
     */
    friend std::ostream& operator<<(std::ostream& os, rand48& rng)
    {
	state_type state_;
	rng.save(state_);
	os << state_.x << " " << state_.y << " " << state_.z << " ";
	return os;
    }

    /**
     * restore generator state from text-mode input stream
     */
    friend std::istream& operator>>(std::istream& is, rand48& rng)
    {
	state_type state_;
	is >> state_.x >> state_.y >> state_.z;
	rng.restore(state_);
	return is;
    }

private:
    cuda::config dim_;
    cuda::vector<ushort3> state_;
};

} // namespace ljgpu

#endif /* ! LJGPU_RNG_RAND48_HPP */
