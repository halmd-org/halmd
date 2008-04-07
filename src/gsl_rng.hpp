/* rng/gsl_rng.hpp
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

#ifndef MDSIM_RNG_GSL_HPP
#define MDSIM_RNG_GSL_HPP

#include <gsl/gsl_rng.h>
#include <cmath>
#include "vector2d.hpp"
#include "vector3d.hpp"


namespace rng { namespace gsl
{

template <const gsl_rng_type*& rng_type>
class rng
{
protected:
    gsl_rng *_rng;

public:
    rng()
    {
	// FIXME set GSL error handler
	if (NULL == (_rng = gsl_rng_alloc(rng_type))) {
	    throw std::bad_alloc();
	}
    }

    rng(const rng<rng_type>& r)
    {
	// FIXME set GSL error handler
	if (NULL == (_rng = gsl_rng_clone(r._rng))) {
	    throw std::bad_alloc();
	}
    }

    ~rng()
    {
	gsl_rng_free(_rng);
    }

    rng<rng_type>& operator=(const rng<rng_type>& r)
    {
	gsl_rng_memcpy(_rng, r._rng);
	return *this;
    }

    /**
     * initialize generator with integer seed
     */
    void set(unsigned long int seed)
    {
	gsl_rng_set(_rng, seed);
    }

    /**
     * generate random integer in generator-dependent interval
     */
    unsigned int get()
    {
	return gsl_rng_get(_rng);
    }

    /**
     * generate random uniform number
     */
    double uniform()
    {
	return gsl_rng_uniform(_rng);
    }

    /**
     * generate 2-dimensional random unit vector
     */
    void unit_vector(vector2d<double>& v)
    {
	double s = 2. * M_PI * uniform();
	v.x = std::cos(s);
	v.y = std::sin(s);
    }

    /**
     * generate 3-dimensional random unit vector
     */
    void unit_vector(vector3d<double>& v)
    {
	//
	// The following method requires an average of 8/Pi =~ 2.55
	// uniform random numbers. It is described in
	//
	// G. Marsaglia, Choosing a Point from the Surface of a Sphere,
	// The Annals of Mathematical Statistics, 1972, 43, 645-646
	//
	// http://projecteuclid.org/euclid.aoms/1177692644#
	//

	double s;

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
typedef rng<gsl_rng_knuthran2002> knuthran2002;
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

}} // namespace rng::gsl

#endif /* ! MDSIM_RNG_GSL_HPP */
