/* CUDA kernel types
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

#ifndef MDSIM_GPU_TYPES_H
#define MDSIM_GPU_TYPES_H


/**
 * multi-float container to facilitate coalesced memory access
 */
union __builtin_align__(16) float432
{
    float4 u4;
    float3 u3;
    float2 u2;
};

static __inline__ __host__ __device__ float432 make_float432(float4 v)
{
    float432 t; t.u4 = v; return t;
}

static __inline__ __host__ __device__ float432 make_float432(float3 v)
{
    float432 t; t.u3 = v; return t;
}

static __inline__ __host__ __device__ float432 make_float432(float2 v)
{
    float432 t; t.u2 = v; return t;
}


#endif /* ! MDSIM_GPU_TYPES_H */
