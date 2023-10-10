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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_YUKAWA_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_YUKAWA_HPP

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
 * define Yukawa potential and parameters
 */
template <typename float_type_>
class yukawa
{
public:
    typedef float_type_ float_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;

    yukawa(
        matrix_type const& amplitude
      , matrix_type const& sigma
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
        float_type dr = sqrt(rr);
        float_type dr_sigma = dr / sigma_(a, b);
        float_type en_pot = amplitude_(a, b) / dr * exp(- dr_sigma);
        float_type fval = en_pot / rr * (1 + dr_sigma);

        return std::make_tuple(fval, en_pot);
    }

    matrix_type const& amplitude() const
    {
        return amplitude_;
    }

    matrix_type const& sigma() const
    {
        return sigma_;
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
    /** interaction strength, in MD units of energy x length*/
    matrix_type amplitude_;
    /** screening length (interaction range parameter), in MD length units */
    matrix_type sigma_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_YUKAWA_HPP */
