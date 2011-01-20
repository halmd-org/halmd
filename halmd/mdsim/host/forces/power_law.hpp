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

#ifndef HALMD_MDSIM_HOST_FORCES_POWER_LAW_HPP
#define HALMD_MDSIM_HOST_FORCES_POWER_LAW_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/tuple/tuple.hpp>
#include <lua.hpp>

#include <halmd/mdsim/host/forces/pair_trunc.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/numeric/pow.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

/**
 * A power-law potential @f$ r^{-n} @f$ is often used for
 * repulsive smooth spheres. A big advantage is
 * its scale invariance (in the absence of a cutoff).
 */

template <typename float_type>
class power_law
{
public:
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;

    static void luaopen(lua_State* L);

    static char const* name() { return "power law"; }
    static char const* module_name() { return "power_law"; }

    power_law(
        unsigned int ntype
      , int index
      , boost::array<float, 3> const& cutoff
      , boost::array<float, 3> const& epsilon
      , boost::array<float, 3> const& sigma
    );

    /**
     * Compute potential and its derivative at squared distance 'rr'
     * for particles of type 'a' and 'b'
     *
     * Call index-dependent template implementations
     * for efficiency of fixed_pow() function.
     */
    boost::tuple<float_type, float_type, float_type> operator() (float_type rr, unsigned a, unsigned b)
    {
        switch (index_) {
            case 6:  return impl_<6>(rr, a, b);
            case 12: return impl_<12>(rr, a, b);
            case 24: return impl_<24>(rr, a, b);
            case 48: return impl_<48>(rr, a, b);
            default:
                LOG_WARNING_ONCE("Using non-optimised force routine for index " << index_);
                return impl_<0>(rr, a, b);
        }
    }

    int const& index() const
    {
        return index_;
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
    /** optimise pow() function by providing the index at compile time
     * @param rr squared distance between particles
     * @param a type of first interacting particle
     * @param b type of second interacting particle
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     * and hypervirial @f$ r \partial_r r \partial_r U(r) @f$
     */
    template <int const_index>
    boost::tuple<float_type, float_type, float_type> impl_(float_type rr, unsigned a, unsigned b)
    {
        // choose arbitrary index_ if template parameter index = 0
        float_type rni;
        if (const_index > 0) {
            rni = fixed_pow<const_index>(sigma_(a, b) / std::sqrt(rr));
        }
        else {
            rni = std::pow(sigma_(a, b) / std::sqrt(rr), index_);
        }
        float_type eps_rni = epsilon_(a, b) * rni;
        float_type fval = (const_index > 0 ? const_index : index_) * eps_rni / rr;
        float_type en_pot = eps_rni - en_cut_(a, b);
        float_type hvir = fixed_pow<2>((const_index > 0 ? const_index : index_)) * eps_rni;

        return boost::make_tuple(fval, en_pot, hvir);
    }

    /** power law index */
    int index_;
    /** interaction strength in MD units */
    matrix_type epsilon_;
    /** interaction range in MD units */
    matrix_type sigma_;
    /** cutoff length in MD units */
    matrix_type r_cut_;
    /** cutoff length in units of sigma */
    matrix_type r_cut_sigma_;
    /** square of cutoff length */
    matrix_type rr_cut_;
    /** potential energy at cutoff in MD units */
    matrix_type en_cut_;
};

}}} // namespace mdsim::host::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_POWER_LAW_HPP */
