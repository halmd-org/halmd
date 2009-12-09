/* Parallelized rand48 random number generator for CUDA
 *
 * Copyright © 2007-2009  Peter Colberg
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

#ifndef HALMD_RNG_RAND48_HPP
#define HALMD_RNG_RAND48_HPP

#include <algorithm>
#include <halmd/algorithm/prefix_sum.hpp>
#include <halmd/rng/gpu/rand48.hpp>
#include <iostream>

namespace halmd
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
    rand48()
      : sym_a(gpu::rand48::a), sym_c(gpu::rand48::c),
        sym_state(gpu::rand48::state) {}

    /**
     * initialize random number generator with CUDA execution dimensions
     */
    rand48(cuda::config const& dim)
      : sym_a(gpu::rand48::a), sym_c(gpu::rand48::c),
        sym_state(gpu::rand48::state), dim_(dim), g_state(dim.threads())
    {
        // copy state pointer to device symbol
        cuda::copy(g_state.data(), sym_state);
    }

    /**
     * initialise device symbols of arbitrary module for rand48 usage
     */
    void init_symbols(cuda::symbol<uint48>& a, cuda::symbol<uint48>& c,
                      cuda::symbol<ushort3*>& state)
    {
        uint48 a_, c_;
        cuda::copy(sym_a, a_);
        cuda::copy(sym_c, c_);
        cuda::copy(a_, a);
        cuda::copy(c_, c);
        cuda::copy(g_state.data(), state);
    }

    /**
     * change random number generator CUDA execution dimensions
     */
    void resize(cuda::config const& dim)
    {
        if (g_state.size() > 0) {
            ushort3 x;
            // save generator state using old dimensions
            save(x);
            // set new CUDA execution dimensions
            dim_ = dim;
            // reallocate global device memory for generator state
            g_state.resize(dim_.threads());
            // copy state pointer to device symbol
            cuda::copy(g_state.data(), sym_state);
            // restore generator state using new dimensions
            restore(x);
        }
        else {
            // set new CUDA execution dimensions
            dim_ = dim;
            // reallocate global device memory for generator state
            g_state.resize(dim_.threads());
            // copy state pointer to device symbol
            cuda::copy(g_state.data(), sym_state);
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
        gpu::rand48::set(g_a, g_c, a, c, seed);
        stream.synchronize();

        // copy leapfrog multiplier into constant device memory
        cuda::copy(a, sym_a);
        // copy leapfrog addend into constant device memory
        cuda::copy(c, sym_c);
    }

    /*
     * fill array with uniform random numbers in [0.0, 1.0)
     */
    void uniform(cuda::vector<float>& r, cuda::stream& stream)
    {
        cuda::configure(dim_.grid, dim_.block, stream);
        gpu::rand48::uniform(r, r.size());
    }

    /**
     * fill array with random integers in [0, 2^32-1]
     */
    void get(cuda::vector<uint>& r, cuda::stream& stream)
    {
        cuda::configure(dim_.grid, dim_.block, stream);
        gpu::rand48::get(r, r.size());
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
        gpu::rand48::save(buf_gpu);
        cuda::copy(buf_gpu, buf, stream);
        stream.synchronize();

        mem = buf[0u];
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
        gpu::rand48::restore(g_a, g_c, a, c, mem);
        stream.synchronize();

        // copy leapfrog multiplier into constant device memory
        cuda::copy(a, sym_a);
        // copy leapfrog addend into constant device memory
        cuda::copy(c, sym_c);
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

private:
    cuda::symbol<uint48>& sym_a;
    cuda::symbol<uint48>& sym_c;
    cuda::symbol<ushort3*>& sym_state;
    cuda::config dim_;
    cuda::vector<ushort3> g_state;
};

} // namespace halmd

#endif /* ! HALMD_RNG_RAND48_HPP */
