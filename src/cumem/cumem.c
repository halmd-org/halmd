/*
 * cumem - Get memory info for given CUDA device
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

#include <cuda/cuda.h>
#include <stdio.h>

int main(int argc, char** argv)
{
    CUcontext ctx;
    CUdevice dev;
    CUresult err;
    int d = 0;
    unsigned int free = 0, total = 0;

    if (argc > 1) {
	d = atoi(argv[1]);
    }

    /* create CUDA context for device */
    if (CUDA_SUCCESS != cuInit(0)) {
	fprintf(stderr, "cumem: failed to initialize the CUDA driver API\n", d);
	return 1;
    }
    if (CUDA_SUCCESS != cuDeviceGet(&dev, d)) {
	fprintf(stderr, "cumem: failed to initialize the CUDA driver API\n", d);
	return 1;
    }
    if (CUDA_SUCCESS != (err = cuCtxCreate(&ctx, 0, dev))) {
	fprintf(stderr, "cumem: failed to create CUDA context for device %d: %d\n", d, err);
	return 1;
    }
    /* query memory info */
    if (CUDA_SUCCESS != cuMemGetInfo(&free, &total)) {
	fprintf(stderr, "cumem: failed to query memory info for device %d\n", d);
	return 1;
    }
    printf("cumem: device %d: %d bytes available, %d bytes total\n", d, free, total);
    return 0;
}
