/*
 * nvlock - Exclusively lock an unused NVIDIA device
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

#include <cuda.h>
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/file.h>

#define NVIDIA_DEVICE_FILENAME "/dev/nvidia%d"

static CUresult (*_cuCtxCreate)(CUcontext *, unsigned int, CUdevice) = 0;
static CUresult (*_cuCtxPopCurrent)(CUcontext *) = 0;
static CUresult (*_cuCtxPushCurrent)(CUcontext) = 0;
static int fd = -1;

CUresult cuCtxCreate(CUcontext *pctx, unsigned int flags, CUdevice dev)
{
    void *handle;
    char fn[16];

    // open dynamic library and load real function symbol
    if (!_cuCtxCreate) {
	handle = dlopen("libcuda.so", RTLD_GLOBAL | RTLD_NOW);
	if (!handle) {
	    return CUDA_ERROR_UNKNOWN;
	}
	_cuCtxCreate = dlsym(handle, "cuCtxCreate");
	if (!_cuCtxCreate) {
	    return CUDA_ERROR_UNKNOWN;
	}
	if (dlclose(handle)) {
	    return CUDA_ERROR_UNKNOWN;
	}
    }

    // lock NVIDIA device file with non-blocking request
    snprintf(fn, sizeof(fn), NVIDIA_DEVICE_FILENAME, dev);
    if (-1 != fd) {
	close(fd);
	fd = -1;
    }
    if (-1 == (fd = open(fn, O_RDWR))) {
	return CUDA_ERROR_UNKNOWN;
    }
    if (-1 == flock(fd, LOCK_EX | LOCK_NB)) {
	close(fd);
	fd = -1;
	return CUDA_ERROR_UNKNOWN;
    }

    // create CUDA context
    return _cuCtxCreate(pctx, flags, dev);
}

CUresult cuCtxPopCurrent(CUcontext *pctx)
{
    void *handle;

    // open dynamic library and load real function symbol
    if (!_cuCtxPopCurrent) {
	handle = dlopen("libcuda.so", RTLD_GLOBAL | RTLD_NOW);
	if (!handle) {
	    return CUDA_ERROR_UNKNOWN;
	}
	_cuCtxPopCurrent = dlsym(handle, "cuCtxPopCurrent");
	if (!_cuCtxPopCurrent) {
	    return CUDA_ERROR_UNKNOWN;
	}
	if (dlclose(handle)) {
	    return CUDA_ERROR_UNKNOWN;
	}
    }

    // unlock the device file
    flock(fd, LOCK_UN);

    // pop current context from CUDA context stack
    return _cuCtxPopCurrent(pctx);
}

CUresult cuCtxPushCurrent(CUcontext ctx)
{
    void *handle;
    CUresult ret;
    CUdevice dev;
    char fn[16];

    // open dynamic library and load real function symbol
    if (!_cuCtxPushCurrent) {
	handle = dlopen("libcuda.so", RTLD_GLOBAL | RTLD_NOW);
	if (!handle) {
	    return CUDA_ERROR_UNKNOWN;
	}
	_cuCtxPushCurrent = dlsym(handle, "cuCtxPushCurrent");
	if (!_cuCtxPushCurrent) {
	    return CUDA_ERROR_UNKNOWN;
	}
	if (dlclose(handle)) {
	    return CUDA_ERROR_UNKNOWN;
	}
    }

    // lock the device file with blocking request
    flock(fd, LOCK_EX);

    // push floating context onto CUDA context stack
    return _cuCtxPushCurrent(ctx);
}
