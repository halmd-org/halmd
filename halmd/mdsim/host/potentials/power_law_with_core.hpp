/*
 * Copyright © 2011  Michael Kopp
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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_POWER_LAW_WITH_CORE_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_POWER_LAW_WITH_CORE_HPP

#include <boost/make_shared.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/tuple/tuple.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/numeric/pow.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {

/**
 * A power-law potential @f$ r^{-n} @f$ is often used for
 * repulsive smooth spheres. This potential is slightly
 * modified to allow for solid particles. It diverges for
 * (finite) r_core (parameter): @f$ (r-r_\text{core})^{-n} @f$
 */

template <typename float_type>
class power_law_with_core
{
public:
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;
    typedef boost::numeric::ublas::symmetric_matrix<unsigned int, boost::numeric::ublas::lower> uint_matrix_type;
    typedef logger logger_type;

    static void luaopen(lua_State* L);

    static char const* module_name() { return "power_law_with_core"; }

    power_law_with_core(
        unsigned int ntype
      , boost::array<float, 3> const& cutoff
      , boost::array<float, 3> const& core
      , boost::array<float, 3> const& epsilon
      , boost::array<float, 3> const& sigma
      , boost::array<unsigned int, 3> const& index
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    /**
     * Compute potential and its derivative at squared distance 'rr'
     * for particles of type 'a' and 'b'
     *
     * Call index-dependent template implementations
     * for efficiency of fixed_pow() function.
     *
     * As of GCC 4.4, the implementation of std::pow has been improved and
     * produces more efficient code than our implementation, fixed_pow.
     *
     * Some benchmarks on a Intel(R) Xeon(R) CPU E5620 @ 2.40GHz
     * using the numeric_pow test:
     *
     * GCC 4.3.2: fixed_pow: 5.29017 ns, std::pow: 9.39675 ns
     * GCC 4.4.1: fixed_pow: 5.10952 ns, std::pow: 1.60335 ns
     *
     * the differences in the assembler code are tiny:
     *
     * std::pow (GCC 4.3.2):
.L738:
        mov     %edx, %eax      # n.3872, n.3872
        addl    $1, %edx        #, n.3872
        cvtsi2sdq       %rax, %xmm0     # n.3872, tmp316
        cmpl    $10000000, %edx #, n.3872
        movapd  %xmm0, %xmm1    # tmp316, tmp324
        mulsd   %xmm0, %xmm1    # tmp316, tmp324
        mulsd   %xmm1, %xmm0    # tmp324, tmp316
        mulsd   %xmm0, %xmm0    # tmp316, tmp316
        mulsd   %xmm0, %xmm0    # tmp316, tmp316
        addsd   %xmm0, %xmm2    # tmp316, b_lsm.3847
        jne     .L738   #,
        movsd   %xmm2, 1304(%rsp)       # b_lsm.3847, b
     *
     * std::pow (GCC 4.4.1):
.L687:
        mov     %eax, %edx      # n.3613, n.3613
        cvtsi2sdq       %rdx, %xmm1     # n.3613, tmp258
        movapd  %xmm1, %xmm0    # tmp258, prephitmp.3563
        mulsd   %xmm1, %xmm0    # tmp258, prephitmp.3563
        mulsd   %xmm1, %xmm0    # tmp258, prephitmp.3563
        mulsd   %xmm0, %xmm0    # prephitmp.3563, prephitmp.3563
        mulsd   %xmm0, %xmm0    # prephitmp.3563, prephitmp.3563
.L683:
        addl    $1, %eax        #, n.3613
        addsd   %xmm0, %xmm2    # prephitmp.3563, prephitmp.3612
        cmpl    $10000000, %eax #, n.3613
        jne     .L687   #,
        movq    72(%rsp), %rdi  # %sfp,
        movsd   %xmm2, 1272(%rsp)       # prephitmp.3612, b
     *
     * fixed_pow (GCC 4.4.1):
.L666:
        mov     %eax, %edx      # n, n
        addl    $1, %eax        #, n
        cvtsi2sdq       %rdx, %xmm0     # n, x.3623
        cmpl    $10000000, %eax #, n
        mulsd   %xmm0, %xmm0    # x.3623, x.3623
        mulsd   %xmm0, %xmm0    # x.3623, x.3623
        movapd  %xmm0, %xmm1    # x.3623, tmp210
        mulsd   %xmm0, %xmm1    # x.3623, tmp210
        mulsd   %xmm0, %xmm1    # x.3623, tmp210
        addsd   %xmm1, %xmm2    # tmp210, prephitmp.3610
        jne     .L666   #,
        leaq    1040(%rsp), %rax        #,
        movsd   %xmm2, 1288(%rsp)       # prephitmp.3610, a
        movq    %rax, %rdi      #,
        movq    %rax, 72(%rsp)  #, %sfp
     *
     *
     */
    boost::tuple<float_type, float_type, float_type> operator() (float_type rr, unsigned a, unsigned b)
    {
        switch (index_(a, b)) {
            case 6:  return impl_<6>(rr, a, b);
            case 12: return impl_<12>(rr, a, b);
            case 24: return impl_<24>(rr, a, b);
            case 48: return impl_<48>(rr, a, b);
            default:
                LOG_WARNING_ONCE("Using non-optimised force routine for index " << index_(a, b));
                return impl_<0>(rr, a, b);
        }
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

    matrix_type const& r_core() const
    {
        return r_core_;
    }

    float_type r_core(unsigned a, unsigned b) const
    {
        return r_core_(a, b);
    }

    float_type rr_core(unsigned a, unsigned b) const
    {
        return rr_core_(a, b);
    }

    matrix_type const& r_core_sigma() const
    {
        return r_core_sigma_;
    }

    matrix_type const& epsilon() const
    {
        return epsilon_;
    }

    matrix_type const& sigma() const
    {
        return sigma_;
    }

    uint_matrix_type const& index() const
    {
        return index_;
    }

private:
    /**
     * Optimise pow() function by providing the index at compile time.
     *
     * @param rr squared distance between particles
     * @param a type of first interacting particle
     * @param b type of second interacting particle
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     * and hypervirial @f$ r \partial_r r \partial_r U(r) @f$ with
     * @f{eqnarray*}{
     *   - \frac{U'(r)}{r} &=& n \epsilon \frac{\sigma^n}{r (r-r_\text{core})^{n+1}} = \frac{n}{r(r-r_\text{core})} U(r) \\
     *   U(r) &=& \epsilon \left(\frac{\sigma}{r - r_\text{core}}\right)^n \\
     *   r \partial_r r \partial_r U(r) &=& n U(r) \left( \frac{r (r_\text{core} + n r) }{ (r-r_\text{core})^2 } \right)
     * @f}
     *
     */
    template <int const_index>
    boost::tuple<float_type, float_type, float_type> impl_(float_type rr, unsigned a, unsigned b)
    {
        // choose arbitrary index_ if template parameter index = 0
        unsigned int n = const_index > 0 ? const_index : index_(a, b);
        float_type rrs = rr / sigma2_(a, b);
        // It's not possible to avoid computation of a squareroot as r_core
        // must be substracted from r but only r*r is passed.
        float_type rs = std::sqrt(rrs);
        float_type rms = rs - r_core_sigma_(a, b);
        float_type rmsi = 1 / rms;

        float_type temp = epsilon_(a, b) * ( (const_index > 0) ?  fixed_pow<const_index+2>(rmsi) : std::pow(rmsi, n+2) );
        float_type en_pot = temp * rms*rms - en_cut_(a, b);
        temp *= n;
        float_type fval = temp * rms / (sigma2_(a, b) * rs);
        float_type hvir = temp * rs * (rms + n * rs);

        /* backup algorithm -- probably slower...
        float_type en_pot = epsilon_(a,b) * ( (const_index > 0) ?  fixed_pow<const_index>(rmsi) : std::pow(rmsi, n) );
        en_pot -= en_cut_(a, b);
        float_type fval = en_pot * n / (sigma2_(a, b) * rs * rms);
        float_type hvir = fval * rr * (r_core_sigma_(a, b) + n * rs)/rms;
        */

        return boost::make_tuple(fval, en_pot, hvir);
    }

    /** interaction strength in MD units */
    matrix_type epsilon_;
    /** interaction range in MD units */
    matrix_type sigma_;
    /** power law index */
    uint_matrix_type index_;
    /** square of pair separation in MD units */
    matrix_type sigma2_;
    /** cutoff length in MD units */
    matrix_type r_cut_;
    /** cutoff length in units of sigma */
    matrix_type r_cut_sigma_;
    /** square of cutoff length */
    matrix_type rr_cut_;
    /** core radius in MD units */
    matrix_type r_core_;
    /** core radius in units of sigma */
    matrix_type r_core_sigma_;
    /** square of core radius */
    matrix_type rr_core_;
    /** potential energy at cutoff in MD units */
    matrix_type en_cut_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;
};

} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_POWER_LAW_WITH_CORE_HPP */
