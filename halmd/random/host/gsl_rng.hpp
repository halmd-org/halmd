/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_RANDOM_HOST_GSL_RNG_HPP
#define HALMD_RANDOM_HOST_GSL_RNG_HPP

#include <gsl/gsl_rng.h>
#include <stdexcept>

namespace halmd
{
namespace random { namespace host
{

/**
 * GSL random number generator C++ wrapper
 *
 * This class models the NumberGenerator, UniformRandomNumberGenerator
 * and (partially) PseudoRandomNumberGenerator concepts in the Boost
 * Random Number library, as well as the RandomNumberGenerator concept
 * in the STL library.
 */
template <gsl_rng_type const* const& rng_type>
class gsl_rng_wrapper
{
public:
    // We use double precision floating-point as the result type,
    // as the gsl_rng_get() function returns 32-bit integer values,
    // whereas some generators are capable of generating numbers
    // with higher precision if sampled vith gsl_rng_uniform().
    typedef double result_type;
    static bool const has_fixed_range = true;

    /**
     * create new instance of random number generator
     */
    gsl_rng_wrapper()
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
     * create new instance of random number generator using given seed
     */
    explicit gsl_rng_wrapper(unsigned long int value)
    {
        // do not abort program after GSL error
        gsl_error_handler_t* handler = gsl_set_error_handler_off();

        if (NULL == (rng_ = gsl_rng_alloc(rng_type))) {
            throw std::bad_alloc();
        }

        // restore previous error handler
        gsl_set_error_handler(handler);

        gsl_rng_set(rng_, value);
    }

    /**
     * create new instance as an exact copy of given generator
     */
    explicit gsl_rng_wrapper(gsl_rng_wrapper const& rng)
    {
        // do not abort program after GSL error
        gsl_error_handler_t* handler = gsl_set_error_handler_off();

        if (NULL == (rng_ = gsl_rng_clone(rng.rng_))) {
            throw std::bad_alloc();
        }

        // restore previous error handler
        gsl_set_error_handler(handler);
    }

    /**
     * copy random number generator into pre-existing generator
     */
    gsl_rng_wrapper& operator=(gsl_rng_wrapper const& rng)
    {
        gsl_rng_memcpy(rng_, rng.rng_);
        return *this;
    }

    /**
     * free all memory associated with generator
     */
    ~gsl_rng_wrapper()
    {
        gsl_rng_free(rng_);
    }

    /**
     * generate uniform pseudo-random floating-point value in [0.0, 1.0)
     */
    result_type operator()()
    {
        return gsl_rng_uniform(rng_);
    }

    /**
     * tight lower bound for pseudo-random floating-point values
     */
    result_type min() const
    {
        return 0.0;
    }

    /**
     * smallest representable number larger than the tight upper bound
     */
    result_type max() const
    {
        return 1.0;
    }

    /**
     * generate uniform integer in [0, n-1]
     */
    unsigned long int operator()(unsigned long int n)
    {
        return gsl_rng_uniform_int(rng_, n);
    }

    /**
     * initialize generator with integer seed
     */
    void seed(unsigned long int value)
    {
        gsl_rng_set(rng_, value);
    }

private:
    gsl_rng *rng_;
};

typedef gsl_rng_wrapper<gsl_rng_borosh13> borosh13;
typedef gsl_rng_wrapper<gsl_rng_coveyou> coveyou;
typedef gsl_rng_wrapper<gsl_rng_cmrg> cmrg;
typedef gsl_rng_wrapper<gsl_rng_fishman18> fishman18;
typedef gsl_rng_wrapper<gsl_rng_fishman20> fishman20;
typedef gsl_rng_wrapper<gsl_rng_fishman2x> fishman2x;
typedef gsl_rng_wrapper<gsl_rng_gfsr4> gfsr4;
typedef gsl_rng_wrapper<gsl_rng_knuthran> knuthran;
typedef gsl_rng_wrapper<gsl_rng_knuthran2> knuthran2;
typedef gsl_rng_wrapper<gsl_rng_lecuyer21> lecuyer21;
typedef gsl_rng_wrapper<gsl_rng_minstd> minstd;
typedef gsl_rng_wrapper<gsl_rng_mrg> mrg;
typedef gsl_rng_wrapper<gsl_rng_mt19937> mt19937;
typedef gsl_rng_wrapper<gsl_rng_mt19937_1999> mt19937_1999;
typedef gsl_rng_wrapper<gsl_rng_mt19937_1998> mt19937_1998;
typedef gsl_rng_wrapper<gsl_rng_r250> r250;
typedef gsl_rng_wrapper<gsl_rng_ran0> ran0;
typedef gsl_rng_wrapper<gsl_rng_ran1> ran1;
typedef gsl_rng_wrapper<gsl_rng_ran2> ran2;
typedef gsl_rng_wrapper<gsl_rng_ran3> ran3;
typedef gsl_rng_wrapper<gsl_rng_rand> rand;
typedef gsl_rng_wrapper<gsl_rng_rand48> rand48;
typedef gsl_rng_wrapper<gsl_rng_random128_bsd> random128_bsd;
typedef gsl_rng_wrapper<gsl_rng_random128_glibc2> random128_glibc2;
typedef gsl_rng_wrapper<gsl_rng_random128_libc5> random128_libc5;
typedef gsl_rng_wrapper<gsl_rng_random256_bsd> random256_bsd;
typedef gsl_rng_wrapper<gsl_rng_random256_glibc2> random256_glibc2;
typedef gsl_rng_wrapper<gsl_rng_random256_libc5> random256_libc5;
typedef gsl_rng_wrapper<gsl_rng_random32_bsd> random32_bsd;
typedef gsl_rng_wrapper<gsl_rng_random32_glibc2> random32_glibc2;
typedef gsl_rng_wrapper<gsl_rng_random32_libc5> random32_libc5;
typedef gsl_rng_wrapper<gsl_rng_random64_bsd> random64_bsd;
typedef gsl_rng_wrapper<gsl_rng_random64_glibc2> random64_glibc2;
typedef gsl_rng_wrapper<gsl_rng_random64_libc5> random64_libc5;
typedef gsl_rng_wrapper<gsl_rng_random8_bsd> random8_bsd;
typedef gsl_rng_wrapper<gsl_rng_random8_glibc2> random8_glibc2;
typedef gsl_rng_wrapper<gsl_rng_random8_libc5> random8_libc5;
typedef gsl_rng_wrapper<gsl_rng_random_bsd> random_bsd;
typedef gsl_rng_wrapper<gsl_rng_random_glibc2> random_glibc2;
typedef gsl_rng_wrapper<gsl_rng_random_libc5> random_libc5;
typedef gsl_rng_wrapper<gsl_rng_randu> randu;
typedef gsl_rng_wrapper<gsl_rng_ranf> ranf;
typedef gsl_rng_wrapper<gsl_rng_ranlux> ranlux;
typedef gsl_rng_wrapper<gsl_rng_ranlux389> ranlux389;
typedef gsl_rng_wrapper<gsl_rng_ranlxd1> ranlxd1;
typedef gsl_rng_wrapper<gsl_rng_ranlxd2> ranlxd2;
typedef gsl_rng_wrapper<gsl_rng_ranlxs0> ranlxs0;
typedef gsl_rng_wrapper<gsl_rng_ranlxs1> ranlxs1;
typedef gsl_rng_wrapper<gsl_rng_ranlxs2> ranlxs2;
typedef gsl_rng_wrapper<gsl_rng_ranmar> ranmar;
typedef gsl_rng_wrapper<gsl_rng_slatec> slatec;
typedef gsl_rng_wrapper<gsl_rng_taus> taus;
typedef gsl_rng_wrapper<gsl_rng_taus2> taus2;
typedef gsl_rng_wrapper<gsl_rng_taus113> taus113;
typedef gsl_rng_wrapper<gsl_rng_transputer> transputer;
typedef gsl_rng_wrapper<gsl_rng_tt800> tt800;
typedef gsl_rng_wrapper<gsl_rng_uni> uni;
typedef gsl_rng_wrapper<gsl_rng_uni32> uni32;
typedef gsl_rng_wrapper<gsl_rng_vax> vax;
typedef gsl_rng_wrapper<gsl_rng_waterman14> waterman14;
typedef gsl_rng_wrapper<gsl_rng_zuf> zuf;

}} // namespace random::host

} // namespace halmd

#endif /* ! HALMD_RANDOM_HOST_GSL_RNG_HPP */
