/**
 * Double-single floating point functions based on DSFLOAT90
 *
 * From the Mandelbrot CUDA programming example by Mark Granger
 * http://forums.nvidia.com/index.php?showtopic=44847
 */

#ifndef HALMD_MATH_DSFLOAT_CUH
#define HALMD_MATH_DSFLOAT_CUH

#include <boost/type_traits/arithmetic_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/math/gpu/dsfun.cuh>

/**
 * double-single floating point value
 */
struct dsfloat
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

/**
 * addition by assignment
 */
__device__ inline dsfloat& operator+=(dsfloat& v, dsfloat const& w)
{
    __dsadd(v.__hi, v.__lo, v.__hi, v.__lo, w.__hi, w.__lo);
    return v;
}

/**
 * subtraction by assignment
 */
__device__ inline dsfloat& operator-=(dsfloat& v, dsfloat const& w)
{
    __dssub(v.__hi, v.__lo, v.__hi, v.__lo, w.__hi, w.__lo);
    return v;
}

/**
 * multiplication by assignment
 */
__device__ inline dsfloat& operator*=(dsfloat& v, dsfloat const& w)
{
    __dsmul(v.__hi, v.__lo, v.__hi, v.__lo, w.__hi, w.__lo);
    return v;
}

/**
 * division by assignment
 */
__device__ inline dsfloat& operator/=(dsfloat& v, dsfloat const& w)
{
    __dsdiv(v.__hi, v.__lo, v.__hi, v.__lo, w.__hi, w.__lo);
    return v;
}

/**
 * addition
 */
__device__ inline dsfloat operator+(dsfloat v, dsfloat const& w)
{
    v += w;
    return v;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator+(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) + w;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator+(dsfloat const& v, T const& w)
{
    return v + static_cast<dsfloat>(w);
}

/**
 * subtraction
 */
__device__ inline dsfloat operator-(dsfloat v, dsfloat const& w)
{
    v -= w;
    return v;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator-(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) - w;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator-(dsfloat const& v, T const& w)
{
    return v - static_cast<dsfloat>(w);
}

/**
 * multiplication
 */
__device__ inline dsfloat operator*(dsfloat v, dsfloat const& w)
{
    v *= w;
    return v;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator*(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) * w;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator*(dsfloat const& v, T const& w)
{
    return v * static_cast<dsfloat>(w);
}

__device__ inline dsfloat operator/(dsfloat v, dsfloat const& w)
{
    v /= w;
    return v;
}

/**
 * division
 */
template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator/(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) / w;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator/(dsfloat const& v, T const& w)
{
    return v / static_cast<dsfloat>(w);
}

/**
 * square root function
 */
__device__ inline dsfloat sqrt(dsfloat v)
{
    dsfloat w;
    __dssqrt(w.__hi, w.__lo, v.__hi, v.__lo);
    return w;
}

#endif /* __CUDACC__ */

#endif /* ! HALMD_MATH_DSFLOAT_CUH */
