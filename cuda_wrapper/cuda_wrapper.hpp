/* cuda_wrapper/cuda_wrapper.hpp
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
 * CUDA runtime API wrapper classes
 */

#ifndef CUDA_WRAPPER_HPP
#define CUDA_WRAPPER_HPP

#ifndef __CUDACC__
#include <cuda_wrapper/error.hpp>

#include <cuda_wrapper/async.hpp>
#include <cuda_wrapper/device.hpp>
#include <cuda_wrapper/thread.hpp>

#include <cuda_wrapper/device/array.hpp>
#include <cuda_wrapper/host/array.hpp>
#endif /* ! __CUDACC__ */

#include <cuda_wrapper/device/function.hpp>
#include <cuda_wrapper/device/symbol.hpp>

#endif /* ! CUDA_WRAPPER_HPP */
