/* cuda_wrapper/texture.hpp
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

#ifndef CUDA_TEXTURE_HPP
#define CUDA_TEXTURE_HPP

#include <cuda/cuda_runtime.h>
#include <cuda_wrapper/error.hpp>
#include <cuda_wrapper/vector.hpp>

/*
 * CUDA texture management
 */

namespace cuda
{

template<class T>
class texture
{
public:
#ifdef __CUDACC__
    /**
     * type-safe constructor for CUDA host code
     */
    texture(::texture<T, 1, cudaReadModeElementType> const& tex) : tex(tex) {}
#endif

    /**
     * bind CUDA texture to device memory array
     */
    void bind(cuda::vector<T> const& g_array)
    {
	CUDA_CALL(cudaBindTexture(NULL, &tex, g_array.data(), &tex.channelDesc));
    }

    /**
     * unbind CUDA texture
     */
    void unbind()
    {
	CUDA_CALL(cudaUnbindTexture(&tex));
    }

private:
    textureReference const& tex;
};

}

#endif /* ! CUDA_TEXTURE_HPP */
