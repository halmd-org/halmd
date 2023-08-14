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

#ifndef HALMD_MDSIM_HOST_MODIFIED_POTENTIALS_PAIR_LENNARD_JONES_HPP
#define HALMD_MDSIM_HOST_MODIFIED_POTENTIALS_PAIR_LENNARD_JONES_HPP

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
 * Define Mie potential and its parameters.
 *
 * @f[ U(r) = 4 \epsilon \left[ (r/\sigma)^{-m} - (r/\sigma)^{-n}) \right] @f]
 *
 * @f$ m, n @f$ must be even and @f$ m > n @f$.
 */
template <typename float_type_>
class mie
{
public:
    typedef float_type_ float_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;
    typedef boost::numeric::ublas::matrix<unsigned> uint_matrix_type;

    mie(
        matrix_type const& epsilon
      , matrix_type const& sigma
      , uint_matrix_type const& index_m
      , uint_matrix_type const& index_n
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /** compute potential and its derivatives at squared distance 'rr' for particles of type 'a' and 'b'
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$, potential @f$ U(r) @f$
     *
     * @f{eqnarray*}{
     *   - U'(r) / r &=& 4 r^{-2} \epsilon (\sigma/r)^{n} \left[ m (\sigma/r)^{m-n} - n \right] \\
     *   U(r) &=& 4 \epsilon (\sigma/r)^{n} \left[ (\sigma/r)^{m-n} - 1 \right]
     * @f}
     */
    std::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
    {
        float_type sigma2 = sigma2_(a, b);
        unsigned m_2 = index_m_2_(a, b);
        unsigned n_2 = index_n_2_(a, b);
        float_type rri = sigma2 / rr;
        float_type rni = pow(rri, n_2);
        float_type rmni = (m_2 - n_2 == n_2) ? rni : pow(rri, m_2 - n_2);
        float_type eps_rni = epsilon_(a, b) * rni;
        float_type fval = 8 * rri * eps_rni * (m_2 * rmni - n_2) / sigma2;
        float_type en_pot = 4 * eps_rni * (rmni - 1);

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

    uint_matrix_type const& index_m() const
    {
        return index_m_;
    }

    uint_matrix_type const& index_n() const
    {
        return index_n_;
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
    /** power law index of repulsion */
    uint_matrix_type index_m_;
    /** half-value of index of repulsion */
    uint_matrix_type index_m_2_;
    /** power law index of attraction */
    uint_matrix_type index_n_;
    /** half-value of index of attraction */
    uint_matrix_type index_n_2_;
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

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_HPP */
