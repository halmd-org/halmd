/**
 * Double-single floating point functions based on DSFUN90
 *
 * From the Mandelbrot CUDA programming example by Mark Granger
 * http://forums.nvidia.com/index.php?showtopic=44847
 */

#ifndef LJGPU_MATH_DSFUN_CUH
#define LJGPU_MATH_DSFUN_CUH

#ifdef __CUDACC__
# include <cuda_runtime.h>
#endif

/**
 * This function sets the DS number A equal to the double precision floating point number B.
 */
__device__ __host__ inline void __dsdeq(float& a0, float& a1, double b)
{
    a0 = (float)b;
    a1 = (float)(b - a0);
}

/**
 * This function sets the DS number A equal to the single precision floating point number B.
 */
__device__ __host__ inline void __dsfeq(float& a0, float& a1, float b)
{
    a0 = b;
    a1 = 0.0f;
}

#ifdef __CUDACC__

/**
 * This function computes c = a + b.
 */
__device__ inline void __dsadd(float& c0, float& c1, float const a0, float const a1, float const b0, float const b1)
{
    // Compute __dsa + __dsb using Knuth's trick.
    float t1 = a0 + b0;
    float e = t1 - a0;
    float t2 = ((b0 - e) + (a0 - (t1 - e))) + a1 + b1;

    // The result is t1 + t2, after normalization.
    c0 = e = t1 + t2;
    c1 = t2 - (e - t1);
}

/**
 * This function computes c = a - b.
 */
__device__ inline void __dssub(float& c0, float& c1, float const a0, float const a1, float const b0, float const b1)
{
    // Compute __dsa - __dsb using Knuth's trick.
    float t1 = a0 - b0;
    float e = t1 - a0;
    float t2 = ((-b0 - e) + (a0 - (t1 - e))) + a1 - b1;

    // The result is t1 + t2, after normalization.
    c0 = e = t1 + t2;
    c1 = t2 - (e - t1);
}

/**
 * This function multiplies DS numbers A and B to yield the DS product C.
 *
 * Modified double-single mul function by Norbert Juffa, NVIDIA
 * uses __fmul_rn() and __fadd_rn() intrinsics which prevent FMAD merging
 *
 * Based on: Guillaume Da Graça, David Defour. Implementation of Float-Float
 * Operators on Graphics Hardware. RNC'7 pp. 23-32, 2006.
 */
__device__ inline void __dsmul(float& c0, float& c1, float const a0, float const a1, float const b0, float const b1)
{
    // This splits __dsa(1) and __dsb(1) into high-order and low-order words.
    float cona = a0 * 8193.0f;
    float conb = b0 * 8193.0f;
    float sa1 = cona - (cona - a0);
    float sb1 = conb - (conb - b0);
    float sa2 = a0 - sa1;
    float sb2 = b0 - sb1;

    // Multiply a0 * b0 using Dekker's method.
    float c11 = __fmul_rn(a0, b0);
    float c21 = (((sa1 * sb1 - c11) + sa1 * sb2) + sa2 * sb1) + sa2 * sb2;

    // Compute a0 * b1 + a1 * b0 (only high-order word is needed).
    float c2 = a0 * b1 + a1 * b0;

    // Compute (c11, c21) + c2 using Knuth's trick, also adding low-order product.
    float t1 = c11 + c2;
    float e = t1 - c11;
    float t2 = ((c2 - e) + (c11 - (t1 - e))) + c21 + a1 * b1;

    // The result is t1 + t2, after normalization.
    c0 = e = t1 + t2;
    c1 = t2 - (e - t1);
}

/**
 * This divides the DS number DSA by the DS number DSB to yield the DS quotient DSC.
 */
__device__ inline void __dsdiv(float& dsc0, float& dsc1, float const dsa0, float const dsa1, float const dsb0, float const dsb1)
{
    // Compute a DP approximation to the quotient.

    float s1 = dsa0 / dsb0;

    // This splits s1 and dsb0 into high-order and low-order words.

    float const split = 8193;
    float cona = s1 * split;
    float conb = dsb0 * split;
    float a1 = cona - (cona - s1);
    float b1 = conb - (conb - dsb0);
    float a2 = s1 - a1;
    float b2 = dsb0 - b1;

    // Multiply s1 * dsb0 using Dekker's method.

    float c11 = s1 * dsb0;
    float c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

    // Compute s1 * dsb1 (only high-order word is needed).

    float c2 = s1 * dsb1;

    // Compute (c11, c21) + c2 using Knuth's trick.

    float t1 = c11 + c2;
    float e = t1 - c11;
    float t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

    // The result is t1 + t2, after normalization.

    float t12 = t1 + t2;
    float t22 = t2 - (t12 - t1);

    // Compute dsa - (t12, t22) using Knuth's trick.

    float t11 = dsa0 - t12;
    e = t11 - dsa0;
    float t21 = ((-t12 - e) + (dsa0 - (t11 - e))) + dsa1 - t22;

    // Compute high-order word of (t11, t21) and divide by dsb0.

    float s2 = (t11 + t21) / dsb0;

    // The result is s1 + s2, after normalization.

    dsc0 = s1 + s2;
    dsc1 = s2 - (dsc0 - s1);
}

/**
 * This subroutine computes dsc = da x db.
 */
__device__ inline void __dsmulss(float& dsc0, float& dsc1, float const da, float const db)
{
    float const split = 8193;

    // This splits da and db into high-order and low-order words.

    float cona = da * split;
    float conb = db * split;
    float a1 = cona - (cona - da);
    float b1 = conb - (conb - db);
    float a2 = da - a1;
    float b2 = db - b1;

    // Multiply da * db using Dekker's method.

    dsc0 = __fmul_rn(da, db);
    dsc1 = (((a1 * b1 - dsc0) + a1 * b2) + a2 * b1) + a2 * b2;
}

/**
 * This computes the square root of the DS number A and returns the DS result in B.
 */
__device__ inline void __dssqrt(float& dsb0, float& dsb1, float const dsa0, float const dsa1)
{
    // This subroutine employs the following formula (due to Alan Karp):
    //
    //       Sqrt(A) = (A * X) + 0.5 * [A - (A * X)^2] * X  (approx.)
    //
    // where X is a double precision approximation to the reciprocal square root,
    // and where the multiplications A * X and [] * X are performed with only
    // double precision.

    if (dsa0 == 0) {
        dsb0 = 0;
        dsb1 = 0;
        return;
    }

    float t1 = 1.f / sqrtf(dsa0);
    float t2 = dsa0 * t1;
    float s00, s01, s10, s11;
    __dsmulss(s00, s01, t2, t2);
    __dssub(s10, s11, dsa0, dsa1, s00, s01);
    float t3 = 0.5f * s10 * t1;
    s00 = t2;
    s01 = 0;
    s10 = t3;
    s11 = 0;
    __dsadd(dsb0, dsb1, s00, s01, s10, s11);
}

#endif /* __CUDACC__ */

#endif /* ! LJGPU_MATH_DSFUN_CUH */
