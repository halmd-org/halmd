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

protected:
    const cuda::config dim;
    /** generator state on device */
    cuda::vector<ushort3> state;

public:
    rand48(cuda::config dim) : dim(dim), state(dim.threads())
    {
    }

    /**
     * initialize generator with 32-bit integer seed
     */
    void init(unsigned int seed, cuda::stream& stream)
    {
	cuda::vector<uint3> A(1);
	cuda::vector<uint3> C(1);

	gpu::rand48::init.configure(dim, stream);
	gpu::rand48::init(&state, &A, &C, seed);
	stream.synchronize();

	// copy leapfrogging multiplier into constant device memory
	gpu::rand48::a = A;
	// copy leapfrogging addend into constant device memory
	gpu::rand48::c = C;
    }

    /**
     * fill array with uniform random numbers
     */
    void uniform(cuda::vector<float>& v, cuda::stream& stream)
    {
	assert(v.size() == dim.threads());
	gpu::rand48::uniform.configure(dim, stream);
	gpu::rand48::uniform(&state, &v, 1);
    }

    /**
     * fill array with 2-dimensional random unit vectors
     */
    void unit_vector(cuda::vector<float2>& v, cuda::stream& stream)
    {
	assert(v.size() == dim.threads());
	gpu::rand48::unit_vector_2d.configure(dim, stream);
	gpu::rand48::unit_vector_2d(&state, &v, 1);
    }

    /**
     * fill array with n-dimensional random unit vectors
     */
    void unit_vector(cuda::vector<float3>& v, cuda::stream& stream)
    {
	assert(v.size() == dim.threads());
	gpu::rand48::unit_vector_3d.configure(dim, stream);
	gpu::rand48::unit_vector_3d(&state, &v, 1);
    }

    /**
     * save generator state to memory
     */
    void save(state_type& mem)
    {
	cuda::vector<ushort3> dbuffer(1);
	cuda::host::vector<ushort3> hbuffer(1);
	cuda::stream stream;

	gpu::rand48::save.configure(dim, stream);
	gpu::rand48::save(&state, &dbuffer);
	hbuffer.memcpy(dbuffer, stream);
	stream.synchronize();

	mem = hbuffer[0];
    }

    /**
     * restore generator state from memory
     */
    void restore(const state_type& mem)
    {
	cuda::vector<uint3> A(1);
	cuda::vector<uint3> C(1);
	cuda::stream stream;

	gpu::rand48::restore.configure(dim, stream);
	gpu::rand48::restore(&state, &A, &C, mem);
	stream.synchronize();

	// copy leapfrogging multiplier into constant device memory
	gpu::rand48::a = A;
	// copy leapfrogging addend into constant device memory
	gpu::rand48::c = C;
    }

    /**
     * save generator state to text-mode output stream
     */
    friend std::ostream& operator<<(std::ostream& os, rand48& rng)
    {
	state_type state;
	rng.save(state);
	os << state.x << " " << state.y << " " << state.z << " ";
	return os;
    }

    /**
     * restore generator state from text-mode input stream
     */
    friend std::istream& operator>>(std::istream& is, rand48& rng)
    {
	state_type state;
	is >> state.x >> state.y >> state.z;
	rng.restore(state);
	return is;
    }
};

} // namespace mdsim

#endif /* ! MDSIM_RAND48_HPP */
