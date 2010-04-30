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

#ifndef HALMD_MDSIM_HOST_RANDOM_HPP
#define HALMD_MDSIM_HOST_RANDOM_HPP

#include <algorithm>
#include <boost/random.hpp>
#include <iterator>

#include <halmd/mdsim/random.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace host
{

class random
  : public mdsim::random
{
public:
    typedef mdsim::random _Base;
    typedef boost::mt19937 random_generator;

    static void resolve(po::options const& vm);
    static void options(po::options_description& desc);
    random(po::options const& vm);
    virtual ~random() {}
    void seed(unsigned int value);

    template <typename input_iterator>
    void shuffle(input_iterator first, input_iterator last);
    template <typename value_type>
    void normal(value_type& r1, value_type& r2, value_type sigma2);

protected:
    /** pseudo-random number generator */
    random_generator rng_;
};

/**
 * Shuffle sequence in-place
 *
 * D.E. Knuth, Art of Computer Programming, Volume 2:
 * Seminumerical Algorithms, 3rd Edition, 1997,
 * Addison-Wesley, pp. 124-125.
 */
template <typename input_iterator>
void random::shuffle(input_iterator first, input_iterator last)
{
    typedef typename std::iterator_traits<input_iterator>::difference_type difference_type;
    boost::uniform_int<difference_type> variate;
    for (difference_type i = last - first; i > 1; --i) {
        std::iter_swap(first + variate(rng_, i), first + (i - 1));
    }
}

/**
 * Generate two random numbers from normal distribution
 *
 * The Box-Muller transformation for generating random numbers
 * in the normal distribution was originally described in
 *
 *   G.E.P. Box and M.E. Muller, A Note on the Generation of
 *   Random Normal Deviates, The Annals of Mathematical Statistics,
 *   1958, 29, p. 610-611
 *
 * Here, we use instead the faster polar method of the Box-Muller
 * transformation, see
 *
 *   D.E. Knuth, Art of Computer Programming, Volume 2: Seminumerical
 *   Algorithms, 3rd Edition, 1997, Addison-Wesley, p. 122
 */
template <typename value_type>
void random::normal(value_type& x, value_type& y, value_type sigma)
{
    boost::uniform_01<random_generator&> variate(rng_);
    value_type s;
    do {
        x = 2. * variate() - 1.;
        y = 2. * variate() - 1.;
        s = x * x + y * y;
    } while (s >= 1.);

    s = sigma * std::sqrt(-2. * std::log(s) / s);
    x *= s;
    y *= s;
}

}} // namespace mdsim::host

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_RANDOM_HPP */
