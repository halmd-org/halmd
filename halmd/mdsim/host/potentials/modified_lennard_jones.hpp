/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_MODIFIED_POTENTIALS_LENNARD_JONES_HPP
#define HALMD_MDSIM_HOST_MODIFIED_POTENTIALS_LENNARD_JONES_HPP

#include <boost/assign.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/tuple/tuple.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {

/**
 * Define modified Lennard-Jones potential and its parameters.
 *
 * @f[ U(r) = 4 \epsilon \left[ (r/\sigma)^{-m} - (r/\sigma)^{-n}) \right] @f]
 *
 * @f$ m, n @f$ must be even and @f$ m > n @f$.
 */
template <typename float_type>
class modified_lennard_jones
{
public:
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;
    typedef boost::numeric::ublas::matrix<unsigned> uint_matrix_type;
    typedef logger logger_type;

    static char const* module_name() { return "modified_lennard_jones"; }

    static void luaopen(lua_State* L);

    modified_lennard_jones(
        unsigned ntype1
      , unsigned ntype2
      , matrix_type const& cutoff
      , matrix_type const& epsilon
      , matrix_type const& sigma
      , uint_matrix_type const& index_m
      , uint_matrix_type const& index_n
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

    /** compute potential and its derivatives at squared distance 'rr' for particles of type 'a' and 'b'
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$, potential @f$ U(r) @f$,
     * and hypervirial @f$ r \partial_r r \partial_r U(r) @f$
     *
     * @f{eqnarray*}{
     *   - U'(r) / r &=& 4 r^{-2} \epsilon (\sigma/r)^{n} \left[ m (\sigma/r)^{m-n} - n \right] \\
     *   U(r) &=& 4 \epsilon (\sigma/r)^{n} \left[ (\sigma/r)^{m-n} - 1 \right] \\
     *   r \partial_r r \partial_r U(r) &=& 4 \epsilon (\sigma/r)^{n} \left[ m^2 (\sigma/r)^{m-n} - n^2 \right]
     * @f}
     */
    boost::tuple<float_type, float_type, float_type> operator() (float_type rr, unsigned a, unsigned b)
    {
        float_type sigma2 = sigma2_(a, b);
        unsigned m_2 = index_m_2_(a, b);
        unsigned n_2 = index_n_2_(a, b);
        float_type rri = sigma2 / rr;
        float_type rni = pow(rri, n_2);
        float_type rmni = (m_2 - n_2 == n_2) ? rni : pow(rri, m_2 - n_2);
        float_type eps_rni = epsilon_(a, b) * rni;
        float_type fval = 8 * rri * eps_rni * (m_2 * rmni - n_2) / sigma2;
        float_type en_pot = 4 * eps_rni * (rmni - 1) - en_cut_(a, b);
        float_type hvir = 16 * eps_rni * (m_2 * m_2 * rmni - n_2 * n_2);

        return boost::make_tuple(fval, en_pot, hvir);
    }

    matrix_type const& r_cut() const
    {
        return r_cut_;
    }

    float_type r_cut(unsigned a, unsigned b) const
    {
        return r_cut_(a, b);
    }

    float_type rr_cut(unsigned a, unsigned b) const
    {
        return rr_cut_(a, b);
    }

    matrix_type const& r_cut_sigma() const
    {
        return r_cut_sigma_;
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
    /** cutoff length in units of sigma */
    matrix_type r_cut_sigma_;
    /** cutoff length in MD units */
    matrix_type r_cut_;
    /** square of cutoff length */
    matrix_type rr_cut_;
    /** square of pair separation */
    matrix_type sigma2_;
    /** potential energy at cutoff length in MD units */
    matrix_type en_cut_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;
};

} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_MODIFIED_LENNARD_JONES_HPP */
