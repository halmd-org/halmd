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
#include <stdio.h>
#include <stdlib.h>
#include <sys/file.h>
#include <unistd.h>

#define ERROR(...) \
    do { \
	fprintf(stderr, "nvlock: " __VA_ARGS__); \
	exit(1); \
    } while(0)

struct nvlock
{
    nvlock()
    {
	int n, i, fd;
	char fn[16];
	CUcontext ctx;
	CUdevice dev;

	if (CUDA_SUCCESS != cuInit(0)) {
	    ERROR("failed to initialize the CUDA driver API\n");
	}
	if (CUDA_SUCCESS != cuDeviceGetCount(&n)) {
	    ERROR("failed to get CUDA device count\n");
	}
	for (i = 0; i < n; ++i) {
	    snprintf(fn, sizeof(fn), "/dev/nvidia%d", i);
	    if (-1 == (fd = open(fn, O_RDWR))) {
		ERROR("could not open device: %s\n", fn);
	    }
	    if (-1 == flock(fd, LOCK_EX | LOCK_NB)) {
		close(fd);
		continue;
	    }
	    if (CUDA_SUCCESS != cuDeviceGet(&dev, i)) {
		ERROR("failed to get CUDA device: %d\n", i);
	    }
	    if (CUDA_SUCCESS != cuCtxCreate(&ctx, 0, dev)) {
		ERROR("failed to create CUDA context for device %d\n", i);
	    }
	    return;
	}
	ERROR("no unused CUDA device found\n");
    }
};

/* singleton object */
static nvlock _nvlock;
