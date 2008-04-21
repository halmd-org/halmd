/* Parallelized rand48 random number generator for CUDA
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
 * FIXME ported to CUDA from GSL
 */

/*
 * FIXME description of algorithm
 */

#include "rand48_glue.hpp"
#include "cutil.h"
#include "rand48.h"
using namespace cuda;


namespace rand48
{

/** leapfrogging multiplier */
static __constant__ uint3 a;
/** leapfrogging addend */
static __constant__ uint3 c;


/**
 * calculate leapfrog rand48 multiplier and addend in order to jump
 * ahead N numbers in the random number sequence with a single step
 */
__inline__ __device__ void leapfrog(uint3& A, uint3& C, unsigned int N)
{
    const uint3 a = make_uint3(0xE66D, 0xDEEC, 0x0005);
    const uint3 c = make_uint3(0x000B, 0, 0);

    const uint3 zero = make_uint3(0, 0, 0);
    const uint3 one = make_uint3(1, 0, 0);
    int k;

    //
    // leapfrog multiplier:
    //   A = a^N mod m
    //
    // leapfrog addend:
    //   C = (c * sum(n = 0..(N-1), a^n)) mod m
    //

    A = one;
    C = one;

    for (k = 1; k < N; k++) {
	A = muladd(a, A, zero);
	C = muladd(C, one, A);
    }

    A = muladd(a, A, zero);
    C = muladd(c, C, zero);
}


/**
 * initialize generator with 32-bit integer seed
 */
__global__ void init(ushort3 *y, uint3 *a, uint3 *c, uint seed)
{
    uint3 A, C;
    ushort3 x;

    leapfrog(A, C, GTID + 1);

    if (GTID + 1 == GTDIM) {
	// store leapfrog constants
	*a = A;
	*c = C;
    }

    if (seed == 0) {
	x.x = 0x330E;
	// default seed
	x.y = 0xABCD;
	x.z = 0x1234;
    }
    else {
	x.x = 0x330E;
	x.y = seed & 0xFFFF;
	x.z = (seed >> 16) & 0xFFFF;
    }

    // generate initial states
    y[GTID] = muladd(A, x, C);
}


/**
 * save generator state
 */
__global__ void save(ushort3 *x, ushort3 *state)
{
    if (GTID == 0) {
	*state = x[0];
    }
}


/**
 * restore generate state
 */
__global__ void restore(ushort3 *x, uint3 *a, uint3 *c, ushort3 state)
{
    uint3 A, C;

    leapfrog(A, C, GTID + 1);

    if (GTID + 1 == GTDIM) {
	// store leapfrog constants
	*a = A;
	*c = C;

	x[0] = state;
    }
    else {
	// generate initial states
	x[GTID + 1] = muladd(A, state, C);
    }
}

/**
 * fill array with uniform random numbers
 */
__global__ void uniform(ushort3* state, float* v, unsigned int count)
{
    ushort3 x = state[GTID];

    for (unsigned int k = 0; k < count; k++) {
	v[GTID + k * GTDIM] = uniform(x);
    }

    state[GTID] = x;
}

} // namespace rand48


namespace mdsim { namespace gpu { namespace rand48
{

symbol<uint3> a(::rand48::a);
symbol<uint3> c(::rand48::c);

function<void (ushort3*, uint3*, uint3*, unsigned int)> init(::rand48::init);
function<void (ushort3*, ushort3*)> save(::rand48::save);
function<void (ushort3*, uint3*, uint3*, ushort3)> restore(::rand48::restore);
function<void (ushort3*, float*, unsigned int)> uniform(::rand48::uniform);

}}} // namespace mdsim::gpu::rand48
