/* cuda_wrapper/texture.hpp
 *
 * Copyright (C) 2007  Peter Colberg
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

namespace cuda
{

template<class T>
class texture
{
public:
#ifndef __CUDACC__
    /**
     * bind CUDA texture to device memory array
     */
    texture(texture const& texture_, vector<T> const& vector_)
      : texref_(texture_.texref_)
    {
        CUDA_CALL(cudaBindTexture(NULL, &texref_, vector_.data(), &texref_.channelDesc));
    }

    /**
     * unbind CUDA texture
     */
    ~texture()
    {
        CUDA_CALL(cudaUnbindTexture(&texref_));
    }
#else
    /**
     * store CUDA texture reference
     */
    texture(::texture<T> const& texref)
      : texref_(texref)
    {}
#endif

private:
    ::textureReference const& texref_;
};

}

#endif /* ! CUDA_TEXTURE_HPP */
