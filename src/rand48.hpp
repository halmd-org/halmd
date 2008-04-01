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
template <typename T>
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
    rand48<T>(cuda::config dim) : dim(dim), state(dim.threads())
    {
    }

    /**
     * initialize generator with 32-bit integer seed
     */
    void init(unsigned int seed, cuda::stream& stream)
    {
	cuda::vector<uint3> A(1);
	cuda::vector<uint3> C(1);

	gpu::rand48<T>::init.configure(dim, stream);
	gpu::rand48<T>::init(&state, &A, &C, seed);
	stream.synchronize();

	// copy leapfrogging multiplier into constant device memory
	gpu::rand48<T>::a = A;
	// copy leapfrogging addend into constant device memory
	gpu::rand48<T>::c = C;
    }

    /**
     * fill array with uniform random numbers
     */
    void get_uniform(cuda::vector<float>& v, cuda::stream& stream)
    {
	assert(v.size() == dim.threads());
	gpu::rand48<T>::get_uniform.configure(dim, stream);
	gpu::rand48<T>::get_uniform(&state, &v, 1);
    }

    /**
     * fill array with n-dimensional random unit vectors
     */
    void get_unit_vector(cuda::vector<T>& v, cuda::stream& stream)
    {
	assert(v.size() == dim.threads());
	gpu::rand48<T>::get_unit_vector.configure(dim, stream);
	gpu::rand48<T>::get_unit_vector(&state, &v, 1);
    }

    /**
     * save generator state to memory
     */
    void save(state_type& mem)
    {
	cuda::vector<ushort3> dbuffer(1);
	cuda::host::vector<ushort3> hbuffer(1);
	cuda::stream stream;

	gpu::rand48<T>::save.configure(dim, stream);
	gpu::rand48<T>::save(&state, &dbuffer);
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

	gpu::rand48<T>::restore.configure(dim, stream);
	gpu::rand48<T>::restore(&state, &A, &C, mem);
	stream.synchronize();

	// copy leapfrogging multiplier into constant device memory
	gpu::rand48<T>::a = A;
	// copy leapfrogging addend into constant device memory
	gpu::rand48<T>::c = C;
    }

    /**
     * save generator state to text-mode output stream
     */
    friend std::ostream& operator<<(std::ostream& os, rand48<T>& rng)
    {
	state_type state;
	rng.save(state);
	os << state.x << " " << state.y << " " << state.z;
	return os;
    }

    /**
     * restore generator state from text-mode input stream
     */
    friend std::istream& operator>>(std::istream& is, rand48<T>& rng)
    {
	state_type state;
	is >> state.x >> std::ws >> state.y >> std::ws >> state.z;
	rng.restore(state);
	return is;
    }
};

} // namespace mdsim

#endif /* ! MDSIM_RAND48_HPP */
