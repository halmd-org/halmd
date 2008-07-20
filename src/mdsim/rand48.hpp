/* Parallelized rand48 random number generator for CUDA
 *
 * Copyright (C) 2007  Peter Colberg
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

#ifndef MDSIM_RAND48_HPP
#define MDSIM_RAND48_HPP

#include <algorithm>
#include <iostream>
#include "gpu/rand48_glue.hpp"
#include "gpu/ljfluid_glue.hpp"


namespace mdsim
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
	ushort3 state;
	bool init = (state_.size() > 0) ? true : false;

	if (init) {
	    // save generator state using old dimensions
	    save(state);
	}

	// set new CUDA execution dimensions
	dim_ = dim;
	// reallocate global device memory for generator state
	state_.resize(dim_.threads());

	if (init) {
	    // restore generator state using new dimensions
	    restore(state);
	}
    }

    /**
     * initialize generator with 32-bit integer seed
     */
    void set(unsigned int seed)
    {
	cuda::vector<uint3> a(1), c(1);

	cuda::configure(dim_.grid, dim_.block);
	gpu::rand48::init(state_, a, c, seed);
	cuda::thread::synchronize();

	// copy leapfrogging multiplier into constant device memory
	cuda::copy(a, gpu::rand48::a);
	cuda::copy(a, mdsim::gpu::ljfluid::a);
	// copy leapfrogging addend into constant device memory
	cuda::copy(c, gpu::rand48::c);
	cuda::copy(c, mdsim::gpu::ljfluid::c);
    }

    /*
     * fill array with uniform random numbers between [0.0, 1.0)
     */
    void uniform(cuda::vector<float>& r, cuda::stream& stream)
    {
	assert(r.size() == dim_.threads());
	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::uniform(state_, r, 1);
    }

    /**
     * generate in-place random permutation of a host array
     */
    template <typename T>
    void shuffle(T& array, cuda::stream& stream)
    {
	assert(array.size() <= dim_.threads());
	// generate uniform random numbers between [0.0, 1.0)
	cuda::vector<float> g_r(dim_.threads());
	cuda::host::vector<float> h_r(dim_.threads());
	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::uniform(state_, g_r, 1);
	// copy random numbers from GPU to host
	cuda::copy(g_r, h_r, stream);
	stream.synchronize();

	//
	// D.E. Knuth, Art of Computer Programming, Volume 2:
	// Seminumerical Algorithms, 3rd Edition, 1997,
	// Addison-Wesley, pp. 124–125.
	//
	for (typename T::size_type n = array.size(); n > 1; --n) {
	    std::swap(array[n * h_r[n - 1]], array[n - 1]);
        }
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
	cuda::vector<uint3> a(1), c(1);
	cuda::stream stream;

	cuda::configure(dim_.grid, dim_.block, stream);
	gpu::rand48::restore(state_, a, c, mem);
	stream.synchronize();

	// copy leapfrogging multiplier into constant device memory
	cuda::copy(a, gpu::rand48::a);
	cuda::copy(a, mdsim::gpu::ljfluid::a);
	// copy leapfrogging addend into constant device memory
	cuda::copy(c, gpu::rand48::c);
	cuda::copy(c, mdsim::gpu::ljfluid::c);
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

    /**
     * get pointer to CUDA device memory
     */
    cuda::vector<ushort3>& state()
    {
	return state_;
    }

private:
    cuda::config dim_;
    cuda::vector<ushort3> state_;
};

} // namespace mdsim

#endif /* ! MDSIM_RAND48_HPP */
