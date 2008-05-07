/* Test alternatives for summing over a CUDA device vector
 *
 * Copyright (C) 2008  Peter Colberg
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

#include <cuda_wrapper.hpp>
#include <stdio.h>
#include "sum_gpu.hpp"
#include "timer.hpp"
#include "math.h"


int main(int argc, char** argv)
{
    cuda::device::set(1);

    printf("# N\tblocks\tthreads\tGPU v_avg\tGPU time [ms]\thost v_avg\thost time [ms]\tspeedup\n");

    for (int j = 1; j <= 32768; j *= 2) {
	cuda::config dim(dim3(j), dim3(512));
	cuda::vector<float> v(dim), v_sum(dim.blocks_per_grid());
	cuda::host::vector<float> v_host(v.size()), v_sum_host(v_sum.size());
	cuda::stream stream;

	double v_avg, t1, t2;
	timer::timer timer;

	v_host = M_PI;
	v.memcpy(v_host, stream);
	stream.synchronize();

	// N, blocks, threads
	printf("%d\t%d\t%d\t", j * 512, j, 512);
	fflush(stdout);

	//
	// GPU benchmark
	//

	timer.start();

	mdsim::test::gpu::blocksum.configure(dim, sizeof(float) * dim.threads_per_block(), stream);
	mdsim::test::gpu::blocksum(v.data(), v_sum.data());
	v_sum_host.memcpy(v_sum, stream);
	stream.synchronize();

	v_avg = 0;
	for (int i = 0; i < v_sum_host.size(); i++) {
	    v_avg += v_sum_host[i];
	}
	v_avg /= v_host.size();

	timer.stop();

	t1 = timer.elapsed();
	timer.reset();

	// GPU v_avg, GPU time [ms]
	printf("%.4f\t%.8f\t", t1 * 1.e3, v_avg);
	fflush(stdout);

	//
	// host benchmark
	//

	timer.start();

	v_host.memcpy(v, stream);
	stream.synchronize();

	//
	// Note that this way of calculating the average value is hazardous
	// using single precision floating-point arithmetic in the case of large
	// vector sizes, as small vector component values are added to a growing
	// sum value.
	//

	v_avg = 0;
	for (int i = 0; i < v_host.size(); i++) {
	    v_avg += v_host[i];
	}
	v_avg /= v_host.size();

	timer.stop();

	t2 = timer.elapsed();
	timer.reset();

	// host v_avg, host time [ms]
	printf("%.4f\t%.8f\t", t2 * 1.e3, v_avg);
	// speedup
	printf("%.2f\n", t2 / t1);
    }

    printf("\n\n");

    return EXIT_SUCCESS;
}
