/*
 * Copyright © 2007-2012  Peter Colberg
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

#ifndef CUDA_TEXTURE_HPP
#define CUDA_TEXTURE_HPP

#include <cuda_runtime.h>

#ifndef __CUDACC__
# include <cuda_wrapper/error.hpp>
# include <cuda_wrapper/vector.hpp>
#endif

/*
 * CUDA texture management
 */

namespace cuda {

template <typename T, int dim = 1, cudaTextureReadMode mode = cudaReadModeElementType>
class texture
{
public:
#ifdef __CUDACC__
    /**
     * type-safe constructor for CUDA host code
     */
    texture(::texture<T, dim, mode> const& tex) : ref_(tex), desc_(tex.channelDesc) {}

    /**
     * variant constructor for CUDA host code
     *
     * For variant textures we need to override the channel desciptor.
     */
    texture(::texture<void, dim, mode> const& tex) : ref_(tex), desc_(cudaCreateChannelDesc<T>()) {}
#else /* ! __CUDACC__ */
    /**
     * bind CUDA texture to device memory array
     */
    void bind(cuda::vector<T> const& array) const
    {
        CUDA_CALL(cudaBindTexture(NULL, &ref_, array.data(), &desc_));
    }

    /**
     * unbind CUDA texture
     */
    void unbind() const
    {
        CUDA_CALL(cudaUnbindTexture(&ref_));
    }
#endif /* ! __CUDACC__ */

private:
#ifndef __CUDACC__
    // provide dummy host constructor to avoid GCC 4.4 warning
    texture() : ref_(textureReference()) {}
#endif

    textureReference const& ref_;
    cudaChannelFormatDesc const desc_;
};

} // namespace cuda

#endif /* ! CUDA_TEXTURE_HPP */
