/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_LENNARD_JONES_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_LENNARD_JONES_HPP

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
 * define Lennard-Jones potential and parameters
 */
template <typename float_type_>
class lennard_jones
{
public:
    typedef float_type_ float_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;

    lennard_jones(
        matrix_type const& epsilon
      , matrix_type const& sigma
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /** compute potential and its derivative at squared distance 'rr' for particles of type 'a' and 'b' */
    std::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
    {
        float_type sigma2 = sigma2_(a, b);
        float_type rri = sigma2 / rr;
        float_type r6i = rri * rri * rri;
        float_type eps_r6i = epsilon_(a, b) * r6i;
        float_type fval = 48 * rri * eps_r6i * (r6i - 0.5) / sigma2;
        float_type en_pot = 4 * eps_r6i * (r6i - 1);

        return std::make_tuple(fval, en_pot);
    }

    /** compute second and third derivatives
     */
    boost::tuple<float_type, float_type> derivatives(float_type rr, unsigned a, unsigned b) const
    {
        float_type sigma2 = sigma2_(a, b);
        float_type rri = 1 / rr;
        float_type r4i = rri * rri;
        float_type sigma2_rri = sigma2 * rri;
        float_type sigma6_r6i = sigma2_rri * sigma2_rri * sigma2_rri;
        float_type eps_sigma6_r10i = epsilon_(a, b) * sigma6_r6i * r4i;

        float_type second_der = 672 * eps_sigma6_r10i * (sigma6_r6i - 6.0/21.0);
        float_type third_der = 10752 * eps_sigma6_r10i * rri * (5.0/28.0 - sigma6_r6i);
        return boost::make_tuple(second_der, third_der);
    }

    matrix_type const& epsilon() const
    {
        return epsilon_;
    }

    matrix_type const& sigma() const
    {
        return sigma_;
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
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** potential well depths in MD units */
    matrix_type epsilon_;
    /** pair separation in MD units */
    matrix_type sigma_;
    /** square of pair separation */
    matrix_type sigma2_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_LENNARD_JONES_HPP */
