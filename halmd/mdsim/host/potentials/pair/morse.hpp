/*
 * Copyright © 2008-2013 Felix Höfling
 * Copyright © 2008-2011 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_MORSE_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_MORSE_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <lua.hpp>
#include <tuple>
#include <memory>

#include <halmd/io/logger.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

/**
 * define Morse potential and parameters
 */
template <typename float_type_>
class morse
{
public:
    typedef float_type_ float_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;

    morse(
        matrix_type const& epsilon
      , matrix_type const& sigma
      , matrix_type const& r_min
      , matrix_type const& distortion
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /**
     * Compute force and potential for interaction.
     *
     * @param rr squared distance between particles
     * @param a type of first interacting particle
     * @param b type of second interacting particle
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     */
    std::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
    {
        float_type B = distortion_(a, b);
        float_type r_sigma = sqrt(rr) / sigma_(a, b) / B;
        float_type exp_dr = exp(r_min_sigma_(a, b) / B - r_sigma);
        float_type B2 = B * B;
        float_type eps_exp_dr = epsilon_(a, b) * exp_dr / (2 * B2 - 1);
        float_type fval = 2 * eps_exp_dr * (exp_dr - B2) * r_sigma / rr;
        float_type en_pot = eps_exp_dr * (exp_dr - 2 * B2);

        return std::make_tuple(fval, en_pot);
    }

    matrix_type const& epsilon() const
    {
        return epsilon_;
    }

    matrix_type const& sigma() const
    {
        return sigma_;
    }

    matrix_type const& r_min() const
    {
        return r_min_;
    }

    matrix_type const& r_min_sigma() const
    {
        return r_min_sigma_;
    }

    matrix_type const& distortion() const
    {
        return distortion_;
    }

    unsigned int size1() const
    {
        return epsilon_.size1();
    }

    unsigned int size2() const
    {
        return epsilon_.size2();
    }

    /**
     * Bind module to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** depths of potential well in MD units */
    matrix_type epsilon_;
    /** width of potential well in MD units */
    matrix_type sigma_;
    /** position of potential well in MD units */
    matrix_type r_min_;
    /** position of potential well in units of sigma */
    matrix_type r_min_sigma_;
    /** distortion factor B */
    matrix_type distortion_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_MORSE_HPP */
