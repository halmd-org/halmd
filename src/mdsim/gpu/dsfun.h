/**
 * Double-single floating point functions based on DSFUN90
 *
 * From the Mandelbrot CUDA programming example by Mark Granger
 * http://forums.nvidia.com/index.php?showtopic=44847
 */

#ifndef DSFUN_H
#define DSFUN_H

/**
 * This function sets the DS number A equal to the double precision floating point number B.
 */
__device__ inline void __dsdeq(float &a0, float &a1, double b)
{
    a0 = (float)b;
    a1 = (float)(b - a0);
}

/**
 * This function sets the DS number A equal to the single precision floating point number B.
 */
__device__ inline void __dsfeq(float &a0, float &a1, float b)
{
    a0 = b;
    a1 = 0.0f;
}

/**
 * This function computes c = a + b.
 */
__device__ inline void __dsadd(float &c0, float &c1, const float a0, const float a1, const float b0, const float b1)
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
__device__ inline void __dssub(float &c0, T &c1, const T a0, const T a1, const T b0, const T b1)
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
__device__ inline void __dsmul(float &c0, T &c1, const float a0, const U a1, const V b0, const V b1)
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
 * Double-single floating point vector types
 */
__device__  struct __align__(8) dfloat
{
    float f0, f1;
#ifdef __cplusplus
    dfloat(float const& f0, float const& f1) : f0(f0), f1(f1) {}
    dfloat(float const& f0) : f0(f0), f1(0) {}
#endif
};

__device__ struct __align__(16) dfloat2
{
    float2 f0, f1;
#ifdef __cplusplus
    dfloat2(float2 const& f0, float2 const& f1) : f0(f0), f1(f1) {}
    dfloat2(float2 const& f0) : f0(f0), f1(make_float2(0, 0)) {}
#endif
};

__device__ struct dfloat3
{
    float3 f0, f1;
#ifdef __cplusplus
    dfloat3(float3 const& f0, float3 const& f1) : f0(f0), f1(f1) {}
    dfloat3(float3 const& f0) : f0(f0), f1(make_float3(0, 0, 0)) {}
#endif
};

__device__ struct dfloat4
{
    float4 f0, f1;
#ifdef __cplusplus
    dfloat4(float4 const& f0, float4 const& f1) : f0(f0), f1(f1) {}
    dfloat4(float4 const& f0) : f0(f0), f1(make_float4(0, 0, 0, 0)) {}
#endif
};

/**
 * Double-single floating point vector addition
 */
__device__ dfloat operator+(dfloat v, dfloat const& w)
{
    __dsadd(v.f0, v.f1, v.f0, v.f1, w.f0, w.f1);
    return v;
}

__device__ dfloat2 operator+(dfloat2 v, dfloat2 const& w)
{
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dsadd(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    return v;
}

__device__ dfloat3 operator+(dfloat3 v, dfloat3 const& w)
{
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dsadd(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dsadd(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    return v;
}

__device__ dfloat4 operator+(dfloat4 v, dfloat4 const& w)
{
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dsadd(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dsadd(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    __dsadd(v.f0.w, v.f1.w, v.f0.w, v.f1.w, w.f0.w, w.f1.w);
    return v;
}

/**
 * Double-single floating point vector assignment by addition
 */
__device__ dfloat& operator+=(dfloat& v, dfloat const& w)
{
    __dsadd(v.f0, v.f1, v.f0, v.f1, w.f0, w.f1);
    return v;
}

__device__ dfloat2& operator+=(dfloat2& v, dfloat2 const& w)
{
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dsadd(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    return v;
}

__device__ dfloat3& operator+=(dfloat3& v, dfloat3 const& w)
{
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dsadd(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dsadd(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    return v;
}

__device__ dfloat4& operator+=(dfloat4& v, dfloat4 const& w)
{
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dsadd(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dsadd(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    __dsadd(v.f0.w, v.f1.w, v.f0.w, v.f1.w, w.f0.w, w.f1.w);
    return v;
}

/**
 * Double-single floating point vector subtraction
 */
__device__ dfloat operator-(dfloat v, dfloat const& w)
{
    __dssub(v.f0, v.f1, v.f0, v.f1, w.f0, w.f1);
    return v;
}

__device__ dfloat2 operator-(dfloat2 v, dfloat2 const& w)
{
    __dssub(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dssub(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    return v;
}

__device__ dfloat3 operator-(dfloat3 v, dfloat3 const& w)
{
    __dssub(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dssub(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dssub(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    return v;
}

__device__ dfloat4 operator-(dfloat4 v, dfloat4 const& w)
{
    __dssub(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dssub(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dssub(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    __dssub(v.f0.w, v.f1.w, v.f0.w, v.f1.w, w.f0.w, w.f1.w);
    return v;
}

/**
 * Double-single floating point vector assignment by subtraction
 */
__device__ dfloat& operator-=(dfloat& v, dfloat const& w)
{
    __dssub(v.f0, v.f1, v.f0, v.f1, w.f0, w.f1);
    return v;
}

__device__ dfloat2& operator-=(dfloat2& v, dfloat2 const& w)
{
    __dssub(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dssub(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    return v;
}

__device__ dfloat3& operator-=(dfloat3& v, dfloat3 const& w)
{
    __dssub(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dssub(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dssub(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    return v;
}

__device__ dfloat4& operator-=(dfloat4& v, dfloat4 const& w)
{
    __dssub(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dssub(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dssub(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    __dssub(v.f0.w, v.f1.w, v.f0.w, v.f1.w, w.f0.w, w.f1.w);
    return v;
}

/**
 * Double-single floating point vector scalar multiplication
 */
__device__ dfloat operator*(dfloat v, dfloat const& w)
{
    __dsmul(v.f0, v.f1, v.f0, v.f1, w.f0, w.f1);
    return v;
}

__device__ dfloat2 operator*(dfloat2 v, dfloat const& w)
{
    __dsmul(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0, w.f1);
    __dsmul(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0, w.f1);
    return v;
}

__device__ dfloat3 operator*(dfloat3 v, dfloat const& w)
{
    __dsmul(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0, w.f1);
    __dsmul(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0, w.f1);
    __dsmul(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0, w.f1);
    return v;
}

__device__ dfloat4 operator*(dfloat4 v, dfloat const& w)
{
    __dsmul(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0, w.f1);
    __dsmul(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0, w.f1);
    __dsmul(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0, w.f1);
    __dsmul(v.f0.w, v.f1.w, v.f0.w, v.f1.w, w.f0, w.f1);
    return v;
}

/**
 * Double-single floating point vector assignment by scalar multiplication
 */
__device__ dfloat& operator*=(dfloat& v, dfloat const& w)
{
    __dsmul(v.f0, v.f1, v.f0, v.f1, w.f0, w.f1);
    return v;
}

__device__ dfloat2& operator*=(dfloat2& v, dfloat const& w)
{
    __dsmul(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0, w.f1);
    __dsmul(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0, w.f1);
    return v;
}

__device__ dfloat3& operator*=(dfloat3& v, dfloat const& w)
{
    __dsmul(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0, w.f1);
    __dsmul(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0, w.f1);
    __dsmul(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0, w.f1);
    return v;
}

__device__ dfloat4& operator*=(dfloat4& v, dfloat const& w)
{
    __dsmul(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0, w.f1);
    __dsmul(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0, w.f1);
    __dsmul(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0, w.f1);
    __dsmul(v.f0.w, v.f1.w, v.f0.w, v.f1.w, w.f0, w.f1);
    return v;
}

#endif /* ! DSFUN_H */
