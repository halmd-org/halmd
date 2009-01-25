/* Accumulator with statistical evaluation functions
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_MATH_ACCUM_GPU_HPP
#define LJGPU_MATH_ACCUM_GPU_HPP

#include <cuda_wrapper.hpp>
#include <ljgpu/math/accum.hpp>

namespace ljgpu { namespace gpu
{

/**
 * Accumulator with statistical evaluation functions
 */
template <typename T, typename _Base = ljgpu::accumulator<T> >
class accumulator : public _Base
{
public:
    typedef cuda::symbol<unsigned int*> count_symbol;
    typedef cuda::symbol<T*> value_symbol;
    typedef cuda::vector<unsigned int> count_vector;
    typedef cuda::vector<T> value_vector;
    typedef cuda::host::vector<unsigned int> count_host_vector;
    typedef cuda::host::vector<T> value_host_vector;

    accumulator(unsigned int blocks) :
	g_n(blocks), g_m(blocks), g_v(blocks),
	h_n(blocks), h_m(blocks), h_v(blocks) {}

    /**
     * copy global device memory pointers to device symbols
     */
    void symbols(count_symbol& n, value_symbol& m, value_symbol& v)
    {
	cuda::copy(g_n.data(), n);
	cuda::copy(g_m.data(), m);
	cuda::copy(g_v.data(), v);
    }

    /**
     * accumulate results from GPU
     */
    accumulator<T, _Base>& operator()()
    {
	// FIXME performance bottleneck
	cuda::copy(g_n, h_n);
	cuda::copy(g_m, h_m);
	cuda::copy(g_v, h_v);

	_Base::clear();
	for (size_t i = 0; i < h_n.size(); ++i) {
	    _Base::operator+=(_Base(h_n[i], h_m[i], h_v[i]));
	}
	return *this;
    }

private:
    count_vector g_n;
    value_vector g_m, g_v;
    count_host_vector h_n;
    value_host_vector h_m, h_v;
};

}} // namespace ljgpu::gpu

#endif /* ! LJGPU_MATH_ACCUM_GPU_HPP */
