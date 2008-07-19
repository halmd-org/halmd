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
__device__ inline void dsdeq(float &a0, float &a1, double b)
{
    a0 = (float)b;
    a1 = (float)(b - a0);
}

/**
 * This function sets the DS number A equal to the single precision floating point number B.
 */
__device__ inline void dsfeq(float &a0, float &a1, float b)
{
    a0 = b;
    a1 = 0.0f;
}

/**
 * This function computes c = a + b.
 */
template <typename T>
__device__ inline void dsadd(T &c0, T &c1, const T a0, const T a1, const T b0, const T b1)
{
    // Compute dsa + dsb using Knuth's trick.
    T t1 = a0 + b0;
    T e = t1 - a0;
    T t2 = ((b0 - e) + (a0 - (t1 - e))) + a1 + b1;

    // The result is t1 + t2, after normalization.
    c0 = e = t1 + t2;
    c1 = t2 - (e - t1);
}

/**
 * This function computes c = a - b.
 */
template <typename T>
__device__ inline void dssub(T &c0, T &c1, const T a0, const T a1, const T b0, const T b1)
{
    // Compute dsa - dsb using Knuth's trick.
    T t1 = a0 - b0;
    T e = t1 - a0;
    T t2 = ((-b0 - e) + (a0 - (t1 - e))) + a1 - b1;

    // The result is t1 + t2, after normalization.
    c0 = e = t1 + t2;
    c1 = t2 - (e - t1);
}

/**
 * This function multiplies DS numbers A and B to yield the DS product C.
 */
template <typename T, typename U, typename V>
__device__ inline void dsmul(T &c0, T &c1, const U a0, const U a1, const V b0, const V b1)
{
    // This splits dsa(1) and dsb(1) into high-order and low-order words.
    U cona = a0 * 8193.0f;
    V conb = b0 * 8193.0f;
    U sa1 = cona - (cona - a0);
    V sb1 = conb - (conb - b0);
    U sa2 = a0 - sa1;
    V sb2 = b0 - sb1;

    // Multiply a0 * b0 using Dekker's method.
    T c11 = a0 * b0;
    T c21 = (((sa1 * sb1 - c11) + sa1 * sb2) + sa2 * sb1) + sa2 * sb2;

    // Compute a0 * b1 + a1 * b0 (only high-order word is needed).
    T c2 = a0 * b1 + a1 * b0;

    // Compute (c11, c21) + c2 using Knuth's trick, also adding low-order product.
    T t1 = c11 + c2;
    T e = t1 - c11;
    T t2 = ((c2 - e) + (c11 - (t1 - e))) + c21 + a1 * b1;

    // The result is t1 + t2, after normalization.
    c0 = e = t1 + t2;
    c1 = t2 - (e - t1);
}

#endif /* ! DSFUN_H */
