/**
 * Double-single floating point functions based on DSFUN90
 *
 * From the Mandelbrot CUDA programming example by Mark Granger
 * http://forums.nvidia.com/index.php?showtopic=44847
 */

#ifndef LJGPU_MATH_DSFUN_CUH
#define LJGPU_MATH_DSFUN_CUH

#include <cuda/cuda_runtime.h>

/**
 * This function sets the DS number A equal to the double precision floating point number B.
 */
__device__ __host__ inline void __dsdeq(float &a0, float &a1, double b)
{
    a0 = (float)b;
    a1 = (float)(b - a0);
}

/**
 * This function sets the DS number A equal to the single precision floating point number B.
 */
__device__ __host__ inline void __dsfeq(float &a0, float &a1, float b)
{
    a0 = b;
    a1 = 0.0f;
}

/**
 * This function computes c = a + b.
 */
__device__ __host__ inline void __dsadd(float &c0, float &c1, const float a0, const float a1, const float b0, const float b1)
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
template <typename T>
__device__ __host__ inline void __dssub(float &c0, T &c1, const T a0, const T a1, const T b0, const T b1)
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
 */
template <typename T, typename U, typename V>
__device__ __host__ inline void __dsmul(float &c0, T &c1, const float a0, const U a1, const V b0, const V b1)
{
    // This splits __dsa(1) and __dsb(1) into high-order and low-order words.
    float cona = a0 * 8193.0f;
    float conb = b0 * 8193.0f;
    float sa1 = cona - (cona - a0);
    float sb1 = conb - (conb - b0);
    float sa2 = a0 - sa1;
    float sb2 = b0 - sb1;

    // Multiply a0 * b0 using Dekker's method.
    float c11 = a0 * b0;
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
 * double-single floating point value
 */
__device__ __host__  struct dfloat
{
    float f0, f1;

    __device__ __host__ inline dfloat() {}

    __device__ __host__ inline dfloat(float f0, float f1) : f0(f0), f1(f1) {}

    __device__ __host__ inline dfloat(float f0) : f0(f0), f1(0) {}

    __device__ __host__ inline operator float() const
    {
	return f0;
    }

    __device__ __host__ inline operator double() const
    {
	return (double) f0 + (double) f1;
    }
};

__device__ __host__ inline dfloat& operator+=(dfloat& v, dfloat const& w)
{
    __dsadd(v.f0, v.f1, v.f0, v.f1, w.f0, w.f1);
    return v;
}

__device__ __host__ inline dfloat& operator-=(dfloat& v, dfloat const& w)
{
    __dssub(v.f0, v.f1, v.f0, v.f1, w.f0, w.f1);
    return v;
}

__device__ __host__ inline dfloat& operator*=(dfloat& v, dfloat const& w)
{
    __dsmul(v.f0, v.f1, v.f0, v.f1, w.f0, w.f1);
    return v;
}

__device__ __host__ inline dfloat operator+(dfloat v, dfloat const& w)
{
    v += w;
    return v;
}

__device__ __host__ inline dfloat operator-(dfloat v, dfloat const& w)
{
    v -= w;
    return v;
}

__device__ __host__ inline dfloat operator*(dfloat v, dfloat const& w)
{
    v *= w;
    return v;
}

#endif /* ! LJGPU_MATH_DSFUN_CUH */
