/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_FORCES_LJ_HPP
#define HALMD_MDSIM_HOST_FORCES_LJ_HPP

#include <boost/assign.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>
#include <utility>

#include <halmd/mdsim/host/forces/pair_short_ranged.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/options.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

/**
 * define Lennard-Jones potential and parameters
 */
template <typename float_type>
class lj_potential
{
public:
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;

    static char const* name() { return "Lennard-Jones"; }
    static char const* module_name() { return "lennard_jones"; }

    static void options(po::options_description& desc);
    static void luaopen(lua_State* L);

    lj_potential(
        unsigned ntype
      , boost::array<float, 3> const& cutoff
      , boost::array<float, 3> const& epsilon
      , boost::array<float, 3> const& sigma
    );

    /** compute potential and its derivative at squared distance 'rr' for particles of type 'a' and 'b' */
    std::pair<float_type, float_type> operator() (float_type rr, unsigned a, unsigned b)
    {
        float_type sigma2 = sigma2_(a, b);
        float_type rri = sigma2 / rr;
        float_type r6i = rri * rri * rri;
        float_type epsilon = epsilon_(a, b);
        float_type fval = 48 * rri * r6i * (r6i - 0.5) * (epsilon / sigma2);
        float_type en_pot = 4 * epsilon * r6i * (r6i - 1) - en_cut_(a, b);

        return std::make_pair(fval, en_pot);
    }

    static boost::array<float, 3> default_cutoff()
    {
        return boost::assign::list_of(2.5f)(2.5f)(2.5f);
    }

    static boost::array<float, 3> default_epsilon()
    {
        return boost::assign::list_of(1.0f)(1.5f)(0.5f);
    }

    static boost::array<float, 3> default_sigma()
    {
        return boost::assign::list_of(1.0f)(0.8f)(0.88f);
    }

    matrix_type const& r_cut() const { return r_cut_; }

    float_type r_cut(unsigned a, unsigned b) const
    {
        return r_cut_(a, b);
    }

    float_type rr_cut(unsigned a, unsigned b) const
    {
        return rr_cut_(a, b);
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
};

}}} // namespace mdsim::host::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_LJ_HPP */
