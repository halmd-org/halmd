/* GSL random number generator C++ wrapper
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#ifndef HALMD_RNG_GSL_RNG_HPP
#define HALMD_RNG_GSL_RNG_HPP

#include <cmath>
#include <gsl/gsl_rng.h>
#include <vector>

#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>

namespace halmd { namespace gsl
{

/**
 * GSL random number generator
 */
template <const gsl_rng_type* const& rng_type>
class rng
{
public:
    /** random number generator state type */
    typedef std::vector<char> state_type;

public:
    /**
     * create new instance of random number generator
     */
    rng()
    {
        // do not abort program after GSL error
        gsl_error_handler_t* handler = gsl_set_error_handler_off();

        if (NULL == (rng_ = gsl_rng_alloc(rng_type))) {
            throw std::bad_alloc();
        }

        // restore previous error handler
        gsl_set_error_handler(handler);
    }

    /**
     * create new instance as an exact copy of given generator
     */
    rng(rng<rng_type> const& src)
    {
        // do not abort program after GSL error
        gsl_error_handler_t* handler = gsl_set_error_handler_off();

        if (NULL == (rng_ = gsl_rng_clone(src.rng_))) {
            throw std::bad_alloc();
        }

        // restore previous error handler
        gsl_set_error_handler(handler);
    }

    /**
     * copy random number generator into pre-existing generator
     */
    rng<rng_type>& operator=(rng<rng_type> const& src)
    {
        gsl_rng_memcpy(rng_, src.rng_);
        return *this;
    }

    /**
     * free all memory associated with generator
     */
    ~rng()
    {
        gsl_rng_free(rng_);
    }

    /**
     * save generator state
     */
    void save(state_type& state) const
    {
        state.resize(rng_type->size);
        memcpy(state.data(), rng_->state, rng_type->size);
    }

    /**
     * restore generator state
     */
    void restore(state_type const& state)
    {
        assert(state.size() == rng_type->size);
        memcpy(rng_->state, state.data(), rng_type->size);
    }

    /**
     * initialize generator with integer seed
     */
    void set(unsigned long int seed)
    {
        gsl_rng_set(rng_, seed);
    }

    /**
     * generate random integer in algorithm-dependent interval
     */
    unsigned long int get()
    {
        return gsl_rng_get(rng_);
    }

    /**
     * determine minimum random integer
     */
    unsigned long int min() const
    {
        return gsl_rng_min(rng_);
    }

    /**
     * determine maximum random integer
     */
    unsigned long int max() const
    {
        return gsl_rng_max(rng_);
    }

    /**
     * generate random uniform number
     */
    double uniform()
    {
        return gsl_rng_uniform(rng_);
    }

    /**
     * generate 2-dimensional random unit vector
     */
    template <typename T>
    void unit_vector(vector<T, 2>& v)
    {
        T s = 2. * M_PI * uniform();
        v.x = std::cos(s);
        v.y = std::sin(s);
    }

    /**
     * generate 3-dimensional random unit vector
     */
    template <typename T>
    void unit_vector(vector<T, 3>& v)
    {
        //
        // The following method requires an average of 8/Pi =~ 2.55
        // uniform random numbers. It is described in
        //
        // G. Marsaglia, Choosing a Point from the Surface of a Sphere,
        // The Annals of Mathematical Statistics, 1972, 43, p. 645-646
        //
        // http://projecteuclid.org/euclid.aoms/1177692644#
        //

        T s;

        do {
            v.x = 2. * uniform() - 1.;
            v.y = 2. * uniform() - 1.;
            s = v.x * v.x + v.y * v.y;
        } while (s >= 1.);

        v.z = 1. - 2. * s;
        s = 2. * std::sqrt(1. - s);
        v.x *= s;
        v.y *= s;
    }

    /**
     * generate 2 random numbers from Gaussian distribution with given variance
     */
    template <typename T>
    void gaussian(T& r1, T& r2, T const& var)
    {
        //
        // The Box-Muller transformation for generating random numbers
        // in the normal distribution was originally described in
        //
        // G.E.P. Box and M.E. Muller, A Note on the Generation of
        // Random Normal Deviates, The Annals of Mathematical Statistics,
        // 1958, 29, p. 610-611
        //
        // Here, we use instead the faster polar method of the Box-Muller
        // transformation, see
        //
        // D.E. Knuth, Art of Computer Programming, Volume 2: Seminumerical
        // Algorithms, 3rd Edition, 1997, Addison-Wesley, p. 122
        //

        T s;

        do {
            r1 = 2. * uniform() - 1.;
            r2 = 2. * uniform() - 1.;
            s = r1 * r1 + r2 * r2;
        } while (s >= 1.);

        s = std::sqrt(-2. * var * std::log(s) / s);
        r1 *= s;
        r2 *= s;
    }

    template <typename T>
    void gaussian(vector<T, 3>& v, T const& var)
    {
        gaussian(v[0], v[1], var);
        gaussian(v[2], v[0], var);
    }

    template <typename T>
    void gaussian(vector<T, 2>& v, T const& var)
    {
        gaussian(v[0], v[1], var);
    }

    /**
     * in-place array shuffle
     */
    template <typename T>
    void shuffle(T& array)
    {
        //
        // D.E. Knuth, Art of Computer Programming, Volume 2:
        // Seminumerical Algorithms, 3rd Edition, 1997,
        // Addison-Wesley, pp. 124-125.
        //
        for (typename T::size_type i = array.size(); i > 1; --i) {
            typename T::size_type r = static_cast<typename T::size_type>(i * uniform());
            std::swap(array[r], array[i - 1]);
        }
    }

private:
    gsl_rng *rng_;
};


typedef rng<gsl_rng_borosh13> borosh13;
typedef rng<gsl_rng_coveyou> coveyou;
typedef rng<gsl_rng_cmrg> cmrg;
typedef rng<gsl_rng_fishman18> fishman18;
typedef rng<gsl_rng_fishman20> fishman20;
typedef rng<gsl_rng_fishman2x> fishman2x;
typedef rng<gsl_rng_gfsr4> gfsr4;
typedef rng<gsl_rng_knuthran> knuthran;
typedef rng<gsl_rng_knuthran2> knuthran2;
typedef rng<gsl_rng_lecuyer21> lecuyer21;
typedef rng<gsl_rng_minstd> minstd;
typedef rng<gsl_rng_mrg> mrg;
typedef rng<gsl_rng_mt19937> mt19937;
typedef rng<gsl_rng_mt19937_1999> mt19937_1999;
typedef rng<gsl_rng_mt19937_1998> mt19937_1998;
typedef rng<gsl_rng_r250> r250;
typedef rng<gsl_rng_ran0> ran0;
typedef rng<gsl_rng_ran1> ran1;
typedef rng<gsl_rng_ran2> ran2;
typedef rng<gsl_rng_ran3> ran3;
typedef rng<gsl_rng_rand> rand;
typedef rng<gsl_rng_rand48> rand48;
typedef rng<gsl_rng_random128_bsd> random128_bsd;
typedef rng<gsl_rng_random128_glibc2> random128_glibc2;
typedef rng<gsl_rng_random128_libc5> random128_libc5;
typedef rng<gsl_rng_random256_bsd> random256_bsd;
typedef rng<gsl_rng_random256_glibc2> random256_glibc2;
typedef rng<gsl_rng_random256_libc5> random256_libc5;
typedef rng<gsl_rng_random32_bsd> random32_bsd;
typedef rng<gsl_rng_random32_glibc2> random32_glibc2;
typedef rng<gsl_rng_random32_libc5> random32_libc5;
typedef rng<gsl_rng_random64_bsd> random64_bsd;
typedef rng<gsl_rng_random64_glibc2> random64_glibc2;
typedef rng<gsl_rng_random64_libc5> random64_libc5;
typedef rng<gsl_rng_random8_bsd> random8_bsd;
typedef rng<gsl_rng_random8_glibc2> random8_glibc2;
typedef rng<gsl_rng_random8_libc5> random8_libc5;
typedef rng<gsl_rng_random_bsd> random_bsd;
typedef rng<gsl_rng_random_glibc2> random_glibc2;
typedef rng<gsl_rng_random_libc5> random_libc5;
typedef rng<gsl_rng_randu> randu;
typedef rng<gsl_rng_ranf> ranf;
typedef rng<gsl_rng_ranlux> ranlux;
typedef rng<gsl_rng_ranlux389> ranlux389;
typedef rng<gsl_rng_ranlxd1> ranlxd1;
typedef rng<gsl_rng_ranlxd2> ranlxd2;
typedef rng<gsl_rng_ranlxs0> ranlxs0;
typedef rng<gsl_rng_ranlxs1> ranlxs1;
typedef rng<gsl_rng_ranlxs2> ranlxs2;
typedef rng<gsl_rng_ranmar> ranmar;
typedef rng<gsl_rng_slatec> slatec;
typedef rng<gsl_rng_taus> taus;
typedef rng<gsl_rng_taus2> taus2;
typedef rng<gsl_rng_taus113> taus113;
typedef rng<gsl_rng_transputer> transputer;
typedef rng<gsl_rng_tt800> tt800;
typedef rng<gsl_rng_uni> uni;
typedef rng<gsl_rng_uni32> uni32;
typedef rng<gsl_rng_vax> vax;
typedef rng<gsl_rng_waterman14> waterman14;
typedef rng<gsl_rng_zuf> zuf;

}} // namespace halmd::gsl

#endif /* ! HALMD_RNG_GSL_RNG_HPP */
