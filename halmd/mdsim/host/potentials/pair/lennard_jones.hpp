/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_LENNARD_JONES_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_LENNARD_JONES_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/tuple/tuple.hpp>
#include <lua.hpp>
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
template <typename float_type>
class lennard_jones
{
public:
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;

    lennard_jones(
        matrix_type const& cutoff
      , matrix_type const& epsilon
      , matrix_type const& sigma
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /** compute potential and its derivative at squared distance 'rr' for particles of type 'a' and 'b' */
    boost::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
    {
        float_type sigma2 = sigma2_(a, b);
        float_type rri = sigma2 / rr;
        float_type r6i = rri * rri * rri;
        float_type eps_r6i = epsilon_(a, b) * r6i;
        float_type fval = 48 * rri * eps_r6i * (r6i - 0.5) / sigma2;
        float_type en_pot = 4 * eps_r6i * (r6i - 1) - en_cut_(a, b);

        return boost::make_tuple(fval, en_pot);
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
    std::shared_ptr<logger> logger_;
};

} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_LENNARD_JONES_HPP */
