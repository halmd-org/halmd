/**
 * Double-single floating point functions based on DSFLOAT90
 *
 * From the Mandelbrot CUDA programming example by Mark Granger
 * http://forums.nvidia.com/index.php?showtopic=44847
 */

#ifndef LJGPU_MATH_DSFLOAT_CUH
#define LJGPU_MATH_DSFLOAT_CUH

#include <ljgpu/math/gpu/dsfun.cuh>

/**
 * double-single floating point value
 */
__device__ __host__  struct dsfloat
{
    float __hi, __lo;

    __device__ __host__ inline dsfloat() {}

#ifdef __CUDACC__
    __device__ inline dsfloat(float __hi, float __lo) : __hi(__hi), __lo(__lo) {}

    __device__ inline dsfloat(float f)
    {
	__dsfeq(__hi, __lo, f);
    }

    __device__ inline operator float() const
    {
	return __hi;
    }
#else
    __host__ inline dsfloat(double d)
    {
	__dsdeq(__hi, __lo, d);
    }

    __host__ inline operator double() const
    {
	return (double) __hi + (double) __lo;
    }
#endif
};

#ifdef __CUDACC__

__device__ inline dsfloat& operator+=(dsfloat& v, dsfloat const& w)
{
    __dsadd(v.__hi, v.__lo, v.__hi, v.__lo, w.__hi, w.__lo);
    return v;
}

__device__ inline dsfloat& operator-=(dsfloat& v, dsfloat const& w)
{
    __dssub(v.__hi, v.__lo, v.__hi, v.__lo, w.__hi, w.__lo);
    return v;
}

__device__ inline dsfloat& operator*=(dsfloat& v, dsfloat const& w)
{
    __dsmul(v.__hi, v.__lo, v.__hi, v.__lo, w.__hi, w.__lo);
    return v;
}

__device__ inline dsfloat& operator/=(dsfloat& v, dsfloat const& w)
{
    __dsdiv(v.__hi, v.__lo, v.__hi, v.__lo, w.__hi, w.__lo);
    return v;
}

__device__ inline dsfloat operator+(dsfloat v, dsfloat const& w)
{
    v += w;
    return v;
}

__device__ inline dsfloat operator-(dsfloat v, dsfloat const& w)
{
    v -= w;
    return v;
}

__device__ inline dsfloat operator*(dsfloat v, dsfloat const& w)
{
    v *= w;
    return v;
}

__device__ inline dsfloat operator/(dsfloat v, dsfloat const& w)
{
    v /= w;
    return v;
}

__device__ inline dsfloat sqrt(dsfloat v)
{
    dsfloat w;
    __dssqrt(w.__hi, w.__lo, v.__hi, v.__lo);
    return w;
}

#endif /* __CUDACC__ */

#endif /* ! LJGPU_MATH_DSFLOAT_CUH */
