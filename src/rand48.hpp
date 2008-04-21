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

#include "gpu/rand48_glue.hpp"
#include "gpu/ljfluid_glue.hpp"
#include <iostream>


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
    rand48(cuda::config const& dim) : dim_(dim), state_(dim.threads())
    {
    }

    /**
     * initialize generator with 32-bit integer seed
     */
    void set(unsigned int seed)
    {
	cuda::vector<uint3> a(1), c(1);

	gpu::rand48::init.configure(dim_);
	gpu::rand48::init(state_.data(), a.data(), c.data(), seed);
	cuda::thread::synchronize();

	// copy leapfrogging multiplier into constant device memory
	gpu::rand48::a.memcpy(a);
	mdsim::gpu::ljfluid::a.memcpy(a);
	// copy leapfrogging addend into constant device memory
	gpu::rand48::c.memcpy(c);
	mdsim::gpu::ljfluid::c.memcpy(c);
    }

    /*
     * fill array with uniform random numbers
     */
    void uniform(cuda::vector<float>& r, cuda::stream& stream)
    {
	assert(r.size() == dim_.threads());
	gpu::rand48::uniform.configure(dim_, stream);
	gpu::rand48::uniform(state_.data(), r.data(), 1);
    }

    /**
     * save generator state to memory
     */
    void save(state_type& mem)
    {
	cuda::stream stream;
	cuda::vector<ushort3> buf_gpu(1);
	cuda::host::vector<ushort3> buf(1);

	gpu::rand48::save.configure(dim_, stream);
	gpu::rand48::save(state_.data(), buf_gpu.data());
	buf.memcpy(buf_gpu, stream);
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

	gpu::rand48::restore.configure(dim_, stream);
	gpu::rand48::restore(state_.data(), a.data(), c.data(), mem);
	stream.synchronize();

	// copy leapfrogging multiplier into constant device memory
	gpu::rand48::a.memcpy(a);
	mdsim::gpu::ljfluid::a.memcpy(a);
	// copy leapfrogging addend into constant device memory
	gpu::rand48::c.memcpy(c);
	mdsim::gpu::ljfluid::c.memcpy(c);
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
    ushort3* data()
    {
	return state_.data();
    }

private:
    const cuda::config dim_;
    cuda::vector<ushort3> state_;
};

} // namespace mdsim

#endif /* ! MDSIM_RAND48_HPP */
