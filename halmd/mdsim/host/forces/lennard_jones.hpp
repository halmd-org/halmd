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

#ifndef HALMD_MDSIM_HOST_FORCES_LENNARD_JONES_HPP
#define HALMD_MDSIM_HOST_FORCES_LENNARD_JONES_HPP

#include <boost/assign.hpp>
#include <boost/make_shared.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/tuple/tuple.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/forces/pair_trunc.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace forces {

/**
 * define Lennard-Jones potential and parameters
 */
template <typename float_type>
class lennard_jones
{
public:
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;
    typedef logger logger_type;

    static char const* name() { return "Lennard-Jones"; }
    static char const* module_name() { return "lennard_jones"; }

    static void luaopen(lua_State* L);

    lennard_jones(
        unsigned ntype
      , boost::array<float, 3> const& cutoff
      , boost::array<float, 3> const& epsilon
      , boost::array<float, 3> const& sigma
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    /** compute potential and its derivative at squared distance 'rr' for particles of type 'a' and 'b' */
    boost::tuple<float_type, float_type, float_type> operator() (float_type rr, unsigned a, unsigned b)
    {
        float_type sigma2 = sigma2_(a, b);
        float_type rri = sigma2 / rr;
        float_type r6i = rri * rri * rri;
        float_type eps_r6i = epsilon_(a, b) * r6i;
        float_type fval = 48 * rri * eps_r6i * (r6i - 0.5) / sigma2;
        float_type en_pot = 4 * eps_r6i * (r6i - 1) - en_cut_(a, b);
        float_type hvir = 576 * eps_r6i * (r6i - 0.25);

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

private:
    /** potential well depths in MD units */
    matrix_type epsilon_;
    /** pair separation in MD units */
    matrix_type sigma_;
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
    boost::shared_ptr<logger_type> logger_;
};

} // namespace mdsim
} // namespace host
} // namespace forces
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_LENNARD_JONES_HPP */
