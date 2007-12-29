/* cuda_mem.h
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

/*
 * CUDA memory management
 */

#ifndef __CUDA_MEMCPY_H__
#define __CUDA_MEMCPY_H__

#include "cuda_base.h"
#include "cuda_error.h"

#ifndef __CUDACC__

struct cuda_mem
{
  /*
   * allocate linear memory on the device
   */
  template <typename T>
  static void malloc(T **ptr, size_t n)
  {
    CUDA_CALL(cudaMalloc(reinterpret_cast<void **>(ptr), n * sizeof(T)));
  }

  /*
   * allocate page-locked host memory accessible to the device
   */
  template <typename T>
  static void malloc_host(T **ptr, size_t n)
  {
    CUDA_CALL(cudaMallocHost(reinterpret_cast<void **>(ptr), n * sizeof(T)));
  }

  /*
   * free device memory
   */
  template <typename T>
  static void free(T *ptr)
  {
    CUDA_CALL(cudaFree(ptr));
  }

  /*
   * free page-locked host memory
   */
  template <typename T>
  static void free_host(T *ptr)
  {
    CUDA_CALL(cudaFreeHost(ptr));
  }

  /*
   * copy from host memory area to host memory area
   */
  template <typename T>
  static void HtoH(T *dst, const T *src, size_t n)
  {
    CUDA_CALL(cudaMemcpy(dst, src, n * sizeof(T), cudaMemcpyHostToHost));
  }

  /*
   * copy from host memory area to device memory area
   */
  template <typename T>
  static void HtoD(T *dst, const T *src, size_t n)
  {
    CUDA_CALL(cudaMemcpy(dst, src, n * sizeof(T), cudaMemcpyHostToDevice));
  }

  /*
   * copy from device memory area to host memory area
   */
  template <typename T>
  static void DtoH(T *dst, const T *src, size_t n)
  {
    CUDA_CALL(cudaMemcpy(dst, src, n * sizeof(T), cudaMemcpyDeviceToHost));
  }

  /*
   * copy from device memory area to device memory area
   */
  template <typename T>
  static void DtoD(T *dst, const T *src, size_t n)
  {
    CUDA_CALL(cudaMemcpy(dst, src, n * sizeof(T), cudaMemcpyDeviceToDevice));
  }

  /*
   * copy from host memory area to constant device memory area
   */
  template <typename T>
  static void HtoS(const T& symbol, const T *src, size_t n)
  {
    CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<const char *>(&symbol), src, n * sizeof(T), 0, cudaMemcpyHostToDevice));
  }

  /*
   * copy from device memory area to constant device memory area
   */
  template <typename T>
  static void DtoS(const T& symbol, const T *src, size_t n)
  {
    CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<const char *>(&symbol), src, n * sizeof(T), 0, cudaMemcpyDeviceToDevice));
  }

  /*
   * copy from constant device memory area to host memory area
   */
  template <typename T>
  static void StoH(T *dst, const T& symbol, size_t n)
  {
    CUDA_CALL(cudaMemcpyFromSymbol(dst, reinterpret_cast<const char *>(&symbol), n * sizeof(T), 0, cudaMemcpyDeviceToHost));
  }

  /*
   * copy from constant device memory area to device memory area
   */
  template <typename T>
  static void StoD(T *dst, const T& symbol, size_t n)
  {
    CUDA_CALL(cudaMemcpyFromSymbol(dst, reinterpret_cast<const char *>(&symbol), n * sizeof(T), 0, cudaMemcpyDeviceToDevice));
  }
};

#endif /* ! __CUDACC__ */

#endif /* ! __CUDA_MEMCPY_H__ */
