/* CUDA memcpy functions
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

#ifndef CUDA_WRAPPER_COPY_HPP
#define CUDA_WRAPPER_COPY_HPP

#include <assert.h>
#include <boost/array.hpp>
#include <cuda_runtime.h>
#include <cuda_wrapper/host/vector.hpp>
#include <cuda_wrapper/stream.hpp>
#include <cuda_wrapper/symbol.hpp>
#include <cuda_wrapper/vector.hpp>
#include <vector>

namespace cuda
{

/**
 * copy from device memory area to device memory area
 */
template <typename T>
void copy(vector<T> const& src, vector<T>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpy(dst.data(), src.data(), src.size() * sizeof(T), cudaMemcpyDeviceToDevice));
}

template <typename T>
void copy(vector<T> const& src, vector<T>& dst, size_t size)
{
    assert(size <= src.capacity());
    assert(size <= dst.capacity());
    CUDA_CALL(cudaMemcpy(dst.data(), src.data(), size * sizeof(T), cudaMemcpyDeviceToDevice));
}

/**
 * copy from host memory area to device memory area
 */
template <typename T, typename Alloc>
void copy(std::vector<T, Alloc> const& src, vector<T>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpy(dst.data(), src.data(), src.size() * sizeof(T), cudaMemcpyHostToDevice));
}

template <typename T, typename Alloc>
void copy(std::vector<T, Alloc> const& src, vector<T>& dst, size_t size)
{
    assert(size <= src.capacity());
    assert(size <= dst.capacity());
    CUDA_CALL(cudaMemcpy(dst.data(), src.data(), size * sizeof(T), cudaMemcpyHostToDevice));
}

/**
 * copy from device memory area to host memory area
 */
template <typename T, typename Alloc>
void copy(vector<T> const& src, std::vector<T, Alloc>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpy(dst.data(), src.data(), src.size() * sizeof(T), cudaMemcpyDeviceToHost));
}

template <typename T, typename Alloc>
void copy(vector<T> const& src, std::vector<T, Alloc>& dst, size_t size)
{
    assert(size <= src.capacity());
    assert(size <= dst.capacity());
    CUDA_CALL(cudaMemcpy(dst.data(), src.data(), size * sizeof(T), cudaMemcpyDeviceToHost));
}

/**
 * copy from host memory area to host memory area
 */
template <typename T, typename Alloc>
void copy(std::vector<T, Alloc> const& src, std::vector<T, Alloc>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpy(dst.data(), src.data(), src.size() * sizeof(T), cudaMemcpyHostToHost));
}

template <typename T, typename Alloc>
void copy(std::vector<T, Alloc> const& src, std::vector<T, Alloc>& dst, size_t size)
{
    assert(size <= src.capacity());
    assert(size <= dst.capacity());
    CUDA_CALL(cudaMemcpy(dst.data(), src.data(), size * sizeof(T), cudaMemcpyHostToHost));
}

/**
 * copy from device symbol to value
 */
template <typename T>
void copy(symbol<T> const& src, T& dst)
{
    assert(src.size() == 1);
    CUDA_CALL(cudaMemcpyFromSymbol(&dst, reinterpret_cast<char const*>(src.data()), src.size() * sizeof(T), 0, cudaMemcpyDeviceToHost));
}

/**
 * copy from device symbol to host memory area
 */
template <typename T, size_t size>
void copy(symbol<T[]> const& src, boost::array<T, size>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyFromSymbol(dst.data(), reinterpret_cast<char const*>(src.data()), src.size() * sizeof(T), 0, cudaMemcpyDeviceToHost));
}

/**
 * copy from device symbol to host memory area
 */
template <typename T, typename Alloc>
void copy(symbol<T[]> const& src, std::vector<T, Alloc>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyFromSymbol(dst.data(), reinterpret_cast<char const*>(src.data()), src.size() * sizeof(T), 0, cudaMemcpyDeviceToHost));
}

/*
 * copy from device symbol to device memory area
 */
template <typename T>
void copy(symbol<T> const& src, vector<T>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyFromSymbol(dst.data(), reinterpret_cast<char const*>(src.data()), src.size() * sizeof(T), 0, cudaMemcpyDeviceToDevice));
}

/*
 * copy from device symbol to device memory area
 */
template <typename T>
void copy(symbol<T[]> const& src, vector<T>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyFromSymbol(dst.data(), reinterpret_cast<char const*>(src.data()), src.size() * sizeof(T), 0, cudaMemcpyDeviceToDevice));
}

/**
 * copy from value to device symbol
 */
template <typename T>
void copy(T const& src, symbol<T>& dst)
{
    assert(1 == dst.size());
    CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<char const*>(dst.data()), &src, sizeof(T), 0, cudaMemcpyHostToDevice));
}

/**
 * copy from host memory area to device symbol
 */
template <typename T, size_t size>
void copy(boost::array<T, size> const& src, symbol<T[]>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<char const*>(dst.data()), src.data(), src.size() * sizeof(T), 0, cudaMemcpyHostToDevice));
}

/**
 * copy from host memory area to device symbol
 */
template <typename T, typename Alloc>
void copy(std::vector<T, Alloc> const& src, symbol<T[]>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<char const*>(dst.data()), src.data(), src.size() * sizeof(T), 0, cudaMemcpyHostToDevice));
}

/**
 * copy from device memory area to device symbol
 */
template <typename T>
void copy(vector<T> const& src, symbol<T>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<char const*>(dst.data()), src.data(), src.size() * sizeof(T), 0, cudaMemcpyDeviceToDevice));
}

/**
 * copy from device memory area to device symbol
 */
template <typename T>
void copy(vector<T> const& src, symbol<T[]>& dst)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<char const*>(dst.data()), src.data(), src.size() * sizeof(T), 0, cudaMemcpyDeviceToDevice));
}

#if (CUDART_VERSION >= 1010)

/**
 * asynchronous copy from device memory area to device memory area
 */
template <typename T>
void copy(vector<T> const& src, vector<T>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyAsync(dst.data(), src.data(), src.size() * sizeof(T), cudaMemcpyDeviceToDevice, stream.data()));
}

template <typename T>
void copy(vector<T> const& src, vector<T>& dst, size_t size, stream& stream)
{
    assert(size <= src.capacity());
    assert(size <= dst.capacity());
    CUDA_CALL(cudaMemcpyAsync(dst.data(), src.data(), size * sizeof(T), cudaMemcpyDeviceToDevice, stream.data()));
}

/**
 * asynchronous copy from host memory area to device memory area
 */
template <typename T>
void copy(host::vector<T> const& src, vector<T>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyAsync(dst.data(), src.data(), src.size() * sizeof(T), cudaMemcpyHostToDevice, stream.data()));
}

template <typename T>
void copy(host::vector<T> const& src, vector<T>& dst, size_t size, stream& stream)
{
    assert(size <= src.capacity());
    assert(size <= dst.capacity());
    CUDA_CALL(cudaMemcpyAsync(dst.data(), src.data(), size * sizeof(T), cudaMemcpyHostToDevice, stream.data()));
}

/**
 * asynchronous copy from host memory area to host memory area
 */
template <typename T>
void copy(host::vector<T> const& src, host::vector<T>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyAsync(dst.data(), src.data(), src.size() * sizeof(T), cudaMemcpyHostToHost, stream.data()));
}

template <typename T>
void copy(host::vector<T> const& src, host::vector<T>& dst, size_t size, stream& stream)
{
    assert(size <= src.capacity());
    assert(size <= dst.capacity());
    CUDA_CALL(cudaMemcpyAsync(dst.data(), src.data(), size * sizeof(T), cudaMemcpyHostToHost, stream.data()));
}

/**
 * asynchronous copy from device memory area to host memory area
 */
template <typename T>
void copy(vector<T> const& src, host::vector<T>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyAsync(dst.data(), src.data(), src.size() * sizeof(T), cudaMemcpyDeviceToHost, stream.data()));
}

template <typename T>
void copy(vector<T> const& src, host::vector<T>& dst, size_t size, stream& stream)
{
    assert(size <= src.capacity());
    assert(size <= dst.capacity());
    CUDA_CALL(cudaMemcpyAsync(dst.data(), src.data(), size * sizeof(T), cudaMemcpyDeviceToHost, stream.data()));
}

#endif /* CUDART_VERSION >= 1010 */

#if (CUDART_VERSION >= 2000)

/**
 * asynchronous copy from device symbol to host memory area
 */
template <typename T>
void copy(symbol<T> const& src, host::vector<T>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyFromSymbolAsync(dst.data(), reinterpret_cast<char const*>(src.data()), src.size() * sizeof(T), 0, cudaMemcpyDeviceToHost, stream.data()));
}

/**
 * asynchronous copy from device symbol to host memory area
 */
template <typename T>
void copy(symbol<T[]> const& src, host::vector<T>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyFromSymbolAsync(dst.data(), reinterpret_cast<char const*>(src.data()), src.size() * sizeof(T), 0, cudaMemcpyDeviceToHost, stream.data()));
}

/*
 * asynchronous copy from device symbol to device memory area
 */
template <typename T>
void copy(symbol<T> const& src, vector<T>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyFromSymbolAsync(dst.data(), reinterpret_cast<char const*>(src.data()), src.size() * sizeof(T), 0, cudaMemcpyDeviceToDevice, stream.data()));
}

/*
 * asynchronous copy from device symbol to device memory area
 */
template <typename T>
void copy(symbol<T[]> const& src, vector<T>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyFromSymbolAsync(dst.data(), reinterpret_cast<char const*>(src.data()), src.size() * sizeof(T), 0, cudaMemcpyDeviceToDevice, stream.data()));
}

/**
 * asynchronous copy from host memory area to device symbol
 */
template <typename T>
void copy(host::vector<T> const& src, symbol<T>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyToSymbolAsync(reinterpret_cast<char const*>(dst.data()), src.data(), src.size() * sizeof(T), 0, cudaMemcpyHostToDevice, stream.data()));
}

/**
 * asynchronous copy from host memory area to device symbol
 */
template <typename T>
void copy(host::vector<T> const& src, symbol<T[]>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyToSymbolAsync(reinterpret_cast<char const*>(dst.data()), src.data(), src.size() * sizeof(T), 0, cudaMemcpyHostToDevice, stream.data()));
}

/**
 * asynchronous copy from device memory area to device symbol
 */
template <typename T>
void copy(vector<T> const& src, symbol<T>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyToSymbolAsync(reinterpret_cast<char const*>(dst.data()), src.data(), src.size() * sizeof(T), 0, cudaMemcpyDeviceToDevice, stream.data()));
}

/**
 * asynchronous copy from device memory area to device symbol
 */
template <typename T>
void copy(vector<T> const& src, symbol<T[]>& dst, stream& stream)
{
    assert(src.size() == dst.size());
    CUDA_CALL(cudaMemcpyToSymbolAsync(reinterpret_cast<char const*>(dst.data()), src.data(), src.size() * sizeof(T), 0, cudaMemcpyDeviceToDevice, stream.data()));
}

#endif /* CUDART_VERSION >= 2000 */

/**
 * fill device memory array with constant byte value
 */
template <typename T>
void memset(vector<T>& array, int const& value)
{
    CUDA_CALL(cudaMemset(array.data(), value, array.size() * sizeof(T)));
}

template <typename T>
void memset(vector<T>& array, int const& value, size_t size)
{
    assert(size <= array.capacity());
    CUDA_CALL(cudaMemset(array.data(), value, size * sizeof(T)));
}

} // namespace cuda

#endif /* ! CUDA_WRAPPER_COPY_HPP */
