/* Lennard-Jones fluid kernel
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

#include <float.h>
#include "ljfluid_glue.hpp"
#include "cutil.h"
#include "types.h"
#include "vector2d.h"
#include "vector3d.h"
#include "rand48.h"
using namespace cuda;


namespace rand48
{

/** leapfrogging multiplier */
static __constant__ uint3 a;
/** leapfrogging addend */
static __constant__ uint3 c;

}


namespace mdsim
{

/** number of particles */
static __constant__ unsigned int npart;

/** simulation timestemp */
static __constant__ float timestep;
/** periodic box length */
static __constant__ float box;

/** squared cutoff length */
static __constant__ float rr_cut;
/** cutoff energy for Lennard-Jones potential at cutoff length */
static __constant__ float en_cut;

#ifdef USE_CELL
/** number of cells per dimension */
static __constant__ unsigned int ncell;
#endif


#ifdef USE_CELL
/**
 * determine cell index for a particle
 */
template <typename T>
__device__ unsigned int compute_cell(T const& r)
{
    //
    // Mapping the positional coordinates of a particle to its corresponding
    // cell index is the most delicate part of the cell lists update.
    // The safest way is to combine round-towards-zero with a successive
    // integer modulo operation, which comes with a performance penalty.
    //
    // As an efficient alternative, we transform the coordinates to the
    // half-open unit interval [0.0, 1.0) and multiply with the number
    // of cells per dimension afterwards.
    //
    const T cell = (__saturatef(r / box) * (1.f - FLT_EPSILON)) * ncell;
#ifdef DIM_3D
    return (uint(cell.z) * ncell + uint(cell.y)) * ncell + uint(cell.x);
#else
    return uint(cell.y) * ncell + uint(cell.x);
#endif
}
#endif /* USE_CELL */

/**
 * Verlet algorithm for integration of equations of motion
 */
template <typename T>
__device__ void verlet_step(T& r, T& rm, T& v, T const& f)
{
    T t = r;
    // update coordinates
    r = 2. * r - rm + f * (timestep * timestep);
    // update velocity
    v = (r - rm) / (2. * timestep);
    // store previous coordinates
    rm = t;
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <typename T>
__device__ void leapfrog_half_step(T& r, T& R, T& v, T const& f)
{
    // half step velocity
    v += f * (timestep / 2);
    // full step coordinates
    T dr = v * timestep;
    // periodically reduced coordinates
    r += dr;
    r -= floorf(r / box) * box;
    // periodically extended coordinates
    R += dr;
}

/**
 * second leapfrog step of integration of equations of motion
 */
template <typename T>
__device__ void leapfrog_full_step(T& v, T const& f)
{
    // full step velocity
    v += f * (timestep / 2);
}

/**
 * calculate particle force using Lennard-Jones potential
 */
template <typename T>
__device__ void compute_force(T const& r1, T const& r2, T& f, float& en, float& virial)
{
    // particle distance vector
    T r = r1 - r2;
    // enforce periodic boundary conditions
    r -= rintf(r / box) * box;
    // squared particle distance
    float rr = r * r;

    // enforce cutoff length
    if (rr >= rr_cut) return;

    // compute Lennard-Jones force in reduced units
    float rri = 1 / rr;
    float ri6 = rri * rri * rri;
    float fval = 48 * rri * ri6 * (ri6 - 0.5f);

    // add contribution to this particle's force only
    f += fval * r;

    // potential energy contribution from this particle
    en += 2 * ri6 * (ri6 - 1) - en_cut;

    // virial equation sum
    virial += 0.5f * fval * rr;
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <typename T>
__global__ void inteq(T* r, T* R, T* v, T* f)
{
    leapfrog_half_step(r[GTID], R[GTID], v[GTID], f[GTID]);
}

#ifdef USE_CELL
/**
 * compute neighbour cell
 */
template <typename T>
__device__ unsigned int compute_neighbour_cell(T const& offset)
{
#ifdef DIM_3D
    // cell belonging to this execution block
    T cell = make_int3(blockIdx.x % ncell, (blockIdx.x / ncell) % ncell, blockIdx.x / ncell / ncell);
    // neighbour cell of this cell
    T neighbour = make_int3((cell.x + ncell + offset.x) % ncell, (cell.y + ncell + offset.y) % ncell, (cell.z + ncell + offset.z) % ncell);

    return (neighbour.z * ncell + neighbour.y) * ncell + neighbour.x;
#else
    // cell belonging to this execution block
    T cell = make_int2(blockIdx.x % ncell, blockIdx.x / ncell);
    // neighbour cell of this cell
    T neighbour = make_int2((cell.x + ncell + offset.x) % ncell, (cell.y + ncell + offset.y) % ncell);

    return neighbour.y * ncell + neighbour.x;
#endif
}

/**
 * compute forces with particles in a neighbour cell
 */
template <unsigned int block_size, bool same_cell, typename T, typename U>
__device__ void compute_cell_forces(T const* g_cell, int const* g_tag, U const& offset, T const& r, T& f, int const& tag, float& en, float& virial)
{
    __shared__ T s_cell[block_size];
    __shared__ int s_tag[block_size];

    // compute cell index
    unsigned int icell = compute_neighbour_cell(offset);

    // load particles coordinates for cell
    s_cell[threadIdx.x] = g_cell[icell * block_size + threadIdx.x];
    s_tag[threadIdx.x] = g_tag[icell * block_size + threadIdx.x];
    __syncthreads();

    for (unsigned int i = 0; i < block_size; ++i) {
	// skip same particle
	if (same_cell && threadIdx.x == i)
	    continue;

	// skip virtual particles
	if (!IS_REAL_PARTICLE(s_tag[i]))
	    continue;

	compute_force(r, s_cell[i], f, en, virial);
    }
    __syncthreads();
}

/**
 * n-dimensional MD simulation step
 */
template <unsigned int block_size, typename T>
__global__ void mdstep(T const* g_r, T* g_v, T* g_f, int const* g_tag, float* g_en, float* g_virial)
{
    // load particle associated with this thread
    T r = g_r[GTID];
    T v = g_v[GTID];
    int tag = g_tag[GTID];

    // potential energy contribution
    float en = 0.;
    // virial equation sum contribution
    float virial = 0.;

    // Lennard-Jones force calculation
#ifdef DIM_3D
    T f = make_float3(0., 0., 0.);
#else
    T f = make_float2(0., 0.);
#endif

    // calculate forces for this and neighbouring cells
#ifdef DIM_3D
    compute_cell_forces<block_size, true>(g_r, g_tag, make_int3( 0,  0,  0), r, f, tag, en, virial);
    // visit 26 neighbour cells
    for (int x = -1; x <= 1; ++x)
	for (int y = -1; y <= 1; ++y)
	    for (int z = -1; z <= 1; ++z)
		if (x != 0 || y != 0 || z != 0)
		    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(x,  y,  z), r, f, tag, en, virial);
#else
    compute_cell_forces<block_size, true>(g_r, g_tag, make_int2( 0,  0), r, f, tag, en, virial);
    // visit 8 neighbour cells
    for (int x = -1; x <= 1; ++x)
	for (int y = -1; y <= 1; ++y)
	    if (x != 0 || y != 0)
		compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(x, y), r, f, tag, en, virial);
#endif

    // second leapfrog step as part of integration of equations of motion
    leapfrog_full_step(v, f);

    // store particle associated with this thread
    g_v[GTID] = v;
    g_f[GTID] = f;
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

#else /* USE_CELL */

/**
 * MD simulation step
 */
template <typename T>
__global__ void mdstep(T* g_r, T* g_v, T* g_f, float* g_en, float* g_virial)
{
    extern __shared__ T s_r[];

    // load particle associated with this thread
    T r = g_r[GTID];
    T v = g_v[GTID];

    // potential energy contribution
    float en = 0.;
    // virial equation sum contribution
    float virial = 0.;

#ifdef DIM_3D
    T f = make_float3(0., 0., 0.);
#else
    T f = make_float2(0., 0.);
#endif

    // iterate over all blocks
    for (unsigned int k = 0; k < gridDim.x; k++) {
	// load positions of particles within block
	s_r[TID] = g_r[k * blockDim.x + TID];

	__syncthreads();

	// iterate over all particles within block
	for (unsigned int j = 0; j < blockDim.x; j++) {
	    // skip identical particle
	    if (blockIdx.x == k && TID == j)
		continue;
	    // skip placeholder particles
	    if (k * blockDim.x + j >= npart)
		continue;

	    // compute Lennard-Jones force with particle
	    compute_force(r, s_r[j], f, en, virial);
	}

	__syncthreads();
    }

    // second leapfrog step of integration of equations of motion
    leapfrog_full_step(v, f);

    // store particle associated with this thread
    g_v[GTID] = v;
    g_f[GTID] = f;
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

#endif /* USE_CELL */

/**
 * place particles on a face centered cubic (FCC) lattice
 */
template <typename T>
__global__ void lattice(T* g_r)
{
    T r;
#ifdef DIM_3D
    // number of particles along 1 lattice dimension
    const unsigned int n = ceilf(cbrtf(npart / 4.));

    // compose primitive vectors from 1-dimensional index
    r.x = ((GTID >> 2) % n) + ((GTID ^ (GTID >> 1)) & 1) / 2.;
    r.y = ((GTID >> 2) / n % n) + (GTID & 1) / 2.;
    r.z = ((GTID >> 2) / n / n) + (GTID & 2) / 4.;
#else
    // number of particles along 1 lattice dimension
    const unsigned int n = ceilf(sqrtf(npart / 2.));

    // compose primitive vectors from 1-dimensional index
    r.x = ((GTID >> 1) % n) + (GTID & 1) / 2.;
    r.y = ((GTID >> 1) / n) + (GTID & 1) / 2.;
#endif
    g_r[GTID] = r * (box / n);
}

/**
 * generate random n-dimensional Maxwell-Boltzmann distributed velocities
 */
template <typename T>
__global__ void boltzmann(T* g_v, float temp, ushort3* g_rng)
{
    T v;
    ushort3 rng = g_rng[GTID];

    rand48::gaussian(v.x, v.y, temp, rng);
#ifdef DIM_3D
    // Box-Muller transformation strictly generates 2 variates at once
    rand48::gaussian(v.y, v.z, temp, rng);
#endif

    g_rng[GTID] = rng;
    g_v[GTID] = v;
}

#ifdef USE_CELL
/**
 * assign particles to cells
 */
template <unsigned int cell_size, typename T>
__global__ void assign_cells(T const* g_part, T* g_cell, int* g_tag)
{
    __shared__ T s_block[cell_size];
    __shared__ int s_icell[cell_size];

    __shared__ T s_cell[cell_size];
    __shared__ int s_tag[cell_size];
    // number of particles in cell
    unsigned int n = 0;

    // mark all particles in cell as virtual particles
    s_tag[threadIdx.x] = VIRTUAL_PARTICLE;
    __syncthreads();

    for (unsigned int i = 0; i < npart; i += cell_size) {
	// load block of particles from global device memory
	T r = g_part[i + threadIdx.x];
	s_block[threadIdx.x] = r;
	s_icell[threadIdx.x] = compute_cell(r);
	__syncthreads();

	if (threadIdx.x == 0) {
	    for (unsigned int j = 0; j < cell_size && (i + j) < npart; j++) {
		if (s_icell[j] == blockIdx.x) {
		    // store particle in cell
		    s_cell[n] = s_block[j];
		    // store particle number
		    s_tag[n] = REAL_PARTICLE(i + j);
		    // increment particle count in cell
		    ++n;
		}
	    }
	}
	__syncthreads();
    }

    // store cell in global device memory
    g_cell[blockIdx.x * cell_size + threadIdx.x] = s_cell[threadIdx.x];
    g_tag[blockIdx.x * cell_size + threadIdx.x] = s_tag[threadIdx.x];
}

/**
 * examine neighbour cell for particles which moved into this block's cell
 */
template <unsigned int cell_size, typename T, typename U>
__device__ void examine_cell(U const& offset, T const* g_ir, T const* g_iR, T const* g_iv, int const* g_itag, T* s_or, T* s_oR, T* s_ov, int* s_otag, unsigned int& npart)
{
    __shared__ T s_ir[cell_size];
    __shared__ T s_iR[cell_size];
    __shared__ T s_iv[cell_size];
    __shared__ int s_itag[cell_size];
    __shared__ unsigned int s_icell[cell_size];

    // compute cell index
    unsigned int icell = compute_neighbour_cell(offset);

    // load particles in cell from global device memory
    T r = g_ir[icell * cell_size + threadIdx.x];
    s_ir[threadIdx.x] = r;
    s_iR[threadIdx.x] = g_iR[icell * cell_size + threadIdx.x];
    s_iv[threadIdx.x] = g_iv[icell * cell_size + threadIdx.x];
    s_itag[threadIdx.x] = g_itag[icell * cell_size + threadIdx.x];
    // compute new cell
    s_icell[threadIdx.x] = compute_cell(r);
    __syncthreads();

    if (threadIdx.x == 0) {
	for (unsigned int j = 0; j < cell_size; j++) {
	    // skip virtual particles
	    if (!IS_REAL_PARTICLE(s_itag[j])) continue;

	    // if particle belongs to this cell
	    if (s_icell[j] == blockIdx.x && npart < cell_size) {
		// store particle in cell
		s_or[npart] = s_ir[j];
		s_oR[npart] = s_iR[j];
		s_ov[npart] = s_iv[j];
		s_otag[npart] = s_itag[j];
		// increment particle count in cell
		++npart;
	    }
	}
    }
    __syncthreads();
}

/**
 * update cells
 */
template <unsigned int cell_size, typename T>
__global__ void update_cells(T const* g_ir, T const* g_iR, T const* g_iv, int const* g_itag, T* g_or, T* g_oR, T* g_ov, int* g_otag)
{
    __shared__ T s_or[cell_size];
    __shared__ T s_oR[cell_size];
    __shared__ T s_ov[cell_size];
    __shared__ int s_otag[cell_size];
    // number of particles in cell
    unsigned int n = 0;

    // mark all particles in cell as virtual particles
    s_otag[threadIdx.x] = VIRTUAL_PARTICLE;
    __syncthreads();

#ifdef DIM_3D
    examine_cell<cell_size>(make_int3( 0,  0,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    // visit 26 neighbour cells
    examine_cell<cell_size>(make_int3(-1,  0,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0,  0, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0,  0, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
#else
    examine_cell<cell_size>(make_int2( 0,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    // visit 8 neighbour cells
    examine_cell<cell_size>(make_int2(-1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2(+1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2( 0, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2( 0, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2(-1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2(-1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2(+1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2(+1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
#endif

    // store cell in global device memory
    g_or[blockIdx.x * cell_size + threadIdx.x] = s_or[threadIdx.x];
    g_oR[blockIdx.x * cell_size + threadIdx.x] = s_oR[threadIdx.x];
    g_ov[blockIdx.x * cell_size + threadIdx.x] = s_ov[threadIdx.x];
    g_otag[blockIdx.x * cell_size + threadIdx.x] = s_otag[threadIdx.x];
}
#endif  /* USE_CELL */

} // namespace mdsim


namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D
function<void (float3*, float3*, float3*, float3*)> inteq(mdsim::inteq);
#ifdef USE_CELL
function<void (float3 const*, float3*, float3*, int const*, float*, float*)> mdstep(mdsim::mdstep<CELL_SIZE>);
function<void (float3 const*, float3*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE>);
function<void (float3 const*, float3 const*, float3 const*, int const*, float3*, float3*, float3*, int*)> update_cells(mdsim::update_cells<CELL_SIZE>);
#else
function<void (float3*, float3*, float3*, float*, float*)> mdstep(mdsim::mdstep);
#endif
function<void (float3*)> lattice(mdsim::lattice);
function<void (float3*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#else
function<void (float2*, float2*, float2*, float2*)> inteq(mdsim::inteq);
#ifdef USE_CELL
function<void (float2 const*, float2*, float2*, int const*, float*, float*)> mdstep(mdsim::mdstep<CELL_SIZE>);
function<void (float2 const*, float2*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE>);
function<void (float2 const*, float2 const*, float2 const*, int const*, float2*, float2*, float2*, int*)> update_cells(mdsim::update_cells<CELL_SIZE>);
#else
function<void (float2*, float2*, float2*, float*, float*)> mdstep(mdsim::mdstep);
#endif
function<void (float2*)> lattice(mdsim::lattice);
function<void (float2*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#endif

symbol<unsigned int> npart(mdsim::npart);
symbol<float> box(mdsim::box);
symbol<float> timestep(mdsim::timestep);
symbol<float> rr_cut(mdsim::rr_cut);
symbol<float> en_cut(mdsim::en_cut);
#ifdef USE_CELL
symbol<unsigned int> ncell(mdsim::ncell);
#endif

symbol<uint3> a(::rand48::a);
symbol<uint3> c(::rand48::c);

}}} // namespace mdsim::gpu::ljfluid
