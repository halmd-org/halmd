/*
 * cufree - Display amount of free and used memory of all CUDA devices
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

int cufree(int d)
{
    CUcontext ctx;
    CUdevice dev;
    CUresult err;
    unsigned int free = 0, total = 0;

    if (CUDA_SUCCESS != cuDeviceGet(&dev, d)) {
	fprintf(stderr, "cufree: failed to initialize the CUDA driver API\n", d);
	return 1;
    }
    /* create CUDA context for device */
    if (CUDA_SUCCESS != (err = cuCtxCreate(&ctx, 0, dev))) {
	fprintf(stderr, "cufree: failed to create CUDA context for device %d: %d\n", d, err);
	return 1;
    }
    /* query memory info */
    if (CUDA_SUCCESS != cuMemGetInfo(&free, &total)) {
	fprintf(stderr, "cufree: failed to query memory info for device %d\n", d);
	return 1;
    }
    if (CUDA_SUCCESS != (err = cuCtxDetach(ctx))) {
	fprintf(stderr, "cufree: failed to detach CUDA context for device %d: %d\n", d, err);
	return 1;
    }
    printf("Device %d: %8d   %8d   %8d\n", d, total / 1024, (total - free) / 1024, free / 1024);
    return 0;
}

int main(int argc, char** argv)
{
    int count, d, err;

    if (CUDA_SUCCESS != cuInit(0)) {
	fprintf(stderr, "cufree: failed to initialize the CUDA driver API\n", d);
	return 1;
    }
    if (CUDA_SUCCESS != cuDeviceGetCount(&count)) {
	fprintf(stderr, "cufree: failed to initialize the CUDA driver API\n", d);
	return 1;
    }
    printf("             total       used       free\n");
    for (d = 0; d < count; ++d) {
	if (0 != (err = cufree(d))) {
	    return err;
	}
    }
    return 0;
}
