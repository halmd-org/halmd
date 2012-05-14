/*
 * Copyright © 2008-2010 Peter Colberg and Felix Höfling
 * Copyright © 2012 Nicolas Höft
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

#ifndef HALMD_MDSIM_FORCES_TRUNC_LOCAL_R4_HPP
#define HALMD_MDSIM_FORCES_TRUNC_LOCAL_R4_HPP

#include <halmd/config.hpp>

#ifndef __CUDACC__
# include <cmath>
# include <lua.hpp>
#endif

namespace halmd {
namespace mdsim {
namespace forces {
namespace trunc {

/**
 * This class implments the smoothing function in the form
 * @f[ S(r) = S(r)= \frac{((r-r_c)/h)^4}{[h^4 + (r-r_c)^4]} \right] @f]
 * where h is the smoothing parameter given in the constructor.
 */
template <typename float_type>
class local_r4
{
public:
#ifndef __CUDACC__
    static void luaopen(lua_State* L);

    local_r4(float_type h) : rri_smooth_(std::pow(h, -2)) {}
#endif

    /**
     * Calculate the smoothing based on the particles distance r,
     * the cutoff distance r_cut, the scalar force |F(r)|/r f_abs, and
     * the potential U(r).
     * While r and r_cur remain unmodified, f_abs and pot will be altered
     * based on the chosed smoothing algorithm
     */
    HALMD_GPU_ENABLED void operator()(float_type r, float_type r_cut, float_type& f_abs, float_type& pot) const;

private:
    float_type rri_smooth_;
};

template <typename float_type>
HALMD_GPU_ENABLED void local_r4<float_type>::operator()(
    float_type r        // absolute particle distnace
  , float_type r_cut    // cutoff radius
  , float_type& f_abs   // F(r) / |r|
  , float_type& pot     // U(r)
) const
{
    float_type dr = r - r_cut;
    float_type x2 = dr * dr * rri_smooth_;
    float_type x4 = x2 * x2;
    float_type x4i = 1 / (1 + x4);
    // smoothing function
    float_type h0_r = x4 * x4i;
    // first derivative
    float_type h1_r = 4 * dr * rri_smooth_ * x2 * x4i * x4i;
    // apply smoothing function to obtain C¹ force function
    f_abs = h0_r * f_abs - h1_r * (pot / r);
    // apply smoothing function to obtain C² potential function
    pot = h0_r * pot;
}

} // namespace trunc
} // namespace forces
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_FORCES_TRUNC_LOCAL_R4_HPP */
