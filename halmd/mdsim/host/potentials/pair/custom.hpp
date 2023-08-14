/*
 * Copyright © 2023 Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_CUSTOM_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_CUSTOM_HPP

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
 * define custom potential and parameters
 */
template <typename float_type_>
class custom
{
public:
    typedef float_type_ float_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;

    // FIXME rename param[2-3] to sensible identifiers for
    // the parameters of the custom potential
    custom(
        matrix_type const& sigma    // one parameter must be named sigma and is used as unit of length in the truncations
      , matrix_type const& param2
      , matrix_type const& param3
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
        // FIXME
        // put here the actual formulas for the potential energy (en_pot) and
        // the force divided by the pair distance (fval).
        // use float_type, sqrt(rr), param2_(a, b), etc.
        float_type fval = - sigma_(a, b) * param2_(a, b);
        float_type en_pot = param3_(a, b) * rr / 2;

        return std::make_tuple(fval, en_pot);
    }

    matrix_type const& sigma() const
    {
        return sigma_;
    }

    // FIXME rename param2
    matrix_type const& param2() const
    {
        return param2_;
    }

    // FIXME rename param3
    matrix_type const& param3() const
    {
        return param3_;
    }

    unsigned int size1() const
    {
        return sigma_.size1();
    }

    unsigned int size2() const
    {
        return sigma_.size2();
    }

    /**
     * Bind module to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** interaction range parameter, in MD units */
    matrix_type sigma_;
    /** FIXME second potential parameter, in MD units */
    matrix_type param2_;
    /** FIXME third potential parameter, in MD units */
    matrix_type param3_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_CUSTOM_HPP */
