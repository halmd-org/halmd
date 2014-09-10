/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Michael Kopp
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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_POWER_LAW_WITH_CORE_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_POWER_LAW_WITH_CORE_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <lua.hpp>
#include <tuple>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/numeric/pow.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

/**
 * A power-law potential @f$ r^{-n} @f$ is often used for
 * repulsive smooth spheres. This potential is slightly
 * modified to allow for solid particles. It diverges for
 * (finite) r_core (parameter): @f$ (r-r_\mathrm{core})^{-n} @f$
 */

template <typename float_type>
class power_law_with_core
{
public:
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;
    typedef boost::numeric::ublas::matrix<unsigned int> uint_matrix_type;

    power_law_with_core(
        matrix_type const& cutoff
      , matrix_type const& core
      , matrix_type const& epsilon
      , matrix_type const& sigma
      , uint_matrix_type const& index
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
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
    std::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
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
    /**
     * Optimise pow() function by providing the index at compile time.
     *
     * @param rr squared distance between particles
     * @param a type of first interacting particle
     * @param b type of second interacting particle
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     *
     * @f{eqnarray*}{
     *   U(r) &=& \epsilon \left(\frac{r - r_\mathrm{core}}{\sigma}\right)^{-n} \\
     *   - \frac{U'(r)}{r} &=& n \frac{n}{r(r-r_\mathrm{core})} U(r)
     * @f}
     *
     */
    template <int const_index>
    std::tuple<float_type, float_type> impl_(float_type rr, unsigned a, unsigned b) const
    {
        // choose arbitrary index_ if template parameter index = 0
        unsigned int n = const_index > 0 ? const_index : index_(a, b);
        float_type rr_ss = rr / sigma2_(a, b);
        // The computation of the square root can not be avoided
        // as r_core must be substracted from r but only r * r is passed.
        float_type r_s = std::sqrt(rr_ss);
        float_type dri = 1 / (r_s - r_core_sigma_(a, b));
        float_type eps_dri_n = epsilon_(a, b) * ((const_index > 0) ? fixed_pow<const_index>(dri) : halmd::pow(dri, n));

        float_type en_pot = eps_dri_n - en_cut_(a, b);
        float_type n_eps_dri_n_1 = n * dri * eps_dri_n;
        float_type fval = n_eps_dri_n_1 / (sigma2_(a, b) * r_s);

        return std::make_tuple(fval, en_pot);
    }

    /** interaction strength in MD units */
    matrix_type epsilon_;
    /** interaction range in MD units */
    matrix_type sigma_;
    /** power law index */
    uint_matrix_type index_;
    /** square of pair separation in MD units */
    matrix_type sigma2_;
    /** cutoff length in units of sigma */
    matrix_type r_cut_sigma_;
    /** cutoff length in MD units */
    matrix_type r_cut_;
    /** square of cutoff length */
    matrix_type rr_cut_;
    /** core radius in units of sigma */
    matrix_type r_core_sigma_;
    /** potential energy at cutoff in MD units */
    matrix_type en_cut_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_POWER_LAW_WITH_CORE_HPP */
