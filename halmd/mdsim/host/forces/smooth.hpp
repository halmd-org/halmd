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

#ifndef HALMD_MDSIM_HOST_FORCES_SMOOTH_HPP
#define HALMD_MDSIM_HOST_FORCES_SMOOTH_HPP

#include <lua.hpp>

#include <halmd/utility/program_options/program_options.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

/**
 * provide functions to make the potential @f$C^2@f$-smooth
 * at the cutoff
 */
template <int dimension, typename float_type>
class smooth
{
public:
    static void options(po::options_description& desc);

    static void luaopen(lua_State* L);

    smooth(double r_smooth);
    void compute(float_type r, float_type dr, float_type& fval, float_type& pot);

protected:
    /** potential smoothing function scale parameter */
    float_type r_smooth_;
    /** inverse squared smoothing parameter */
    float_type rri_smooth_;
};

/*
 * compute C²-smooth potential
 *
 * FIXME smoothing function from preprint
 */
template <int dimension, typename float_type>
void smooth<dimension, float_type>::compute(
    float_type r        // absolute particle distnace
  , float_type r_cut    // cutoff radius
  , float_type& fval    // F(r) / |r|
  , float_type& en_pot  // U(r)
)
{
    float_type dr = r - r_cut;
    float_type x2 = dr * dr * rri_smooth_;
    float_type x4 = x2 * x2;
    float_type x4i = 1 / (1 + x4);
    // smoothing function
    float_type h0_r = x4 * x4i;
    // first derivative times (r_smooth)^(-1) [sic!]
    float_type h1_r = 4 * dr * rri_smooth_ * x2 * x4i * x4i;
    // apply smoothing function to obtain C¹ force function
    fval = h0_r * fval - h1_r * (en_pot / r);
    // apply smoothing function to obtain C² potential function
    en_pot = h0_r * en_pot;
}

}}} // namespace mdsim::host::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_SMOOTH_HPP */
