/* CUDA device symbol
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

#ifndef CUDA_SYMBOL_HPP
#define CUDA_SYMBOL_HPP

#include <cuda/cuda_runtime.h>
#ifndef __CUDACC__
#include <cuda_wrapper/error.hpp>
#endif

namespace cuda
{

/**
 * CUDA device symbol
 */
template <typename T>
class symbol
{
public:
    typedef symbol<T> vector_type;
    typedef T value_type;
    typedef size_t size_type;

public:
    /**
     * initialize with device symbol variable
     */
    symbol(const value_type& symbol) : size_(1), ptr_(&symbol)
    {
    }

    /**
     * initialize with device symbol vector
     */
    symbol(const value_type* symbol) : size_(0), ptr_(symbol)
    {
    }

#ifndef __CUDACC__

    /**
     * return element count of device symbol vector
     */
    size_type size() const
    {
	if (!size_) {
	    /*
	     * It would be preferable to issue the following CUDA runtime
	     * call directly upon construction. However, the constructor
	     * has to be compilable by the NVIDIA CUDA compiler as well,
	     * which does not support C++ runtime functionality, e.g.
	     * exceptions.
	     */
	    CUDA_CALL(cudaGetSymbolSize(&size_, reinterpret_cast<char const*>(data())));
	    size_ /= sizeof(value_type);
	}

	return size_;
    }

#endif /* ! __CUDACC__ */

    /**
     * returns device pointer to device symbol
     */
    const value_type* data() const
    {
	return ptr_;
    }

private:
    // disable default copy constructor
    symbol(vector_type const&);
    // disable default assignment operator
    vector_type& operator=(vector_type const&);

private:
    /** symbol size */
    mutable size_type size_;
    /** symbol device pointer */
    const value_type* ptr_;
};

} // namespace cuda

#endif /* ! CUDA_SYMBOL_HPP */
