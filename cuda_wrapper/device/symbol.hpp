/* cuda_wrapper/device/symbol.hpp
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

#ifndef CUDA_DEVICE_SYMBOL_HPP
#define CUDA_DEVICE_SYMBOL_HPP

#include <cuda_runtime.h>
#ifndef __CUDACC__
#include <cuda_wrapper/device/array.hpp>
#include <cuda_wrapper/host/array.hpp>
#endif

namespace cuda
{

namespace host
{

#ifndef __CUDACC__
template <typename T>
class array;
#endif

}

namespace device
{

#ifndef __CUDACC__
template <typename T>
class array;
#endif


template <typename T>
class symbol
{
protected:
    const T *ptr;

public:
    symbol(const T& symbol) : ptr(&symbol) {}

#ifndef __CUDACC__
    symbol& operator=(const T& value)
    {
	// copy from host memory area to device symbol
	CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<const char *>(ptr), &value, sizeof(T), 0, cudaMemcpyHostToDevice));
	return *this;
    }

    symbol& operator=(const host::array<T>& array)
    {
	assert(array.dim() == 1);
	// copy from host memory area to device symbol
	CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<const char *>(ptr), array.get_ptr(), sizeof(T), 0, cudaMemcpyHostToDevice));
	return *this;
    }

    symbol& operator=(const array<T>& array)
    {
	assert(array.dim() == 1);
	// copy from device memory area to device symbol
	CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<const char *>(ptr), array.get_ptr(), sizeof(T), 0, cudaMemcpyDeviceToDevice));
	return *this;
    }
#endif /* ! __CUDACC__ */

    size_t dim() const
    {
	return 1;
    }

    T *get_ptr() const
    {
	return ptr;
    }
};

}

}

#endif /* ! CUDA_DEVICE_SYMBOL_HPP */
