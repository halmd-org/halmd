/*
 * Copyright © 2011-2013  Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#ifndef HALMD_OBSERVABLES_SAMPLES_DENSITY_MODE_HPP
#define HALMD_OBSERVABLES_SAMPLES_DENSITY_MODE_HPP

#include <complex>
#include <lua.hpp>

#include <halmd/utility/raw_array.hpp>

namespace halmd {
namespace observables {
namespace samples {

/**
 * Data structure for storing Fourier modes of the particle density.
 *
 * @f$ \rho_{\vec q} = \sum_{i=1}^N \exp(\textrm{i}\vec q \cdot \vec r_i) @f$
 */
class density_mode
{
public:
    typedef raw_array<std::complex<double>> mode_array_type;

    /**
     * Construct sample of given size.
     *
     * @param nq total number of wavevectors
     */
    density_mode(std::size_t nq) : rho_(nq) {}

    /**
     * Returns const reference to density modes, one entry per wavevector.
     *
     * nested list of density modes, @f$ rho[i][j] = \rho_{\vec q}^{(i)} = @f$
     * for wavevector @f$ \vec q = wavevector[j] @f$ and particle types @f$ i @f$
     */
    mode_array_type const& rho() const
    {
        return rho_;
    }

    /**
     * Returns non-const reference to density modus, one entry per wavevector.
     *
     * nested list of density modes, @f$ rho[i][j] = \rho_{\vec q}^{(i)} = @f$
     * for wavevector @f$ \vec q = wavevector[j] @f$ and particle types @f$ i @f$
     */
    mode_array_type& rho()
    {
        return rho_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** density modes */
    mode_array_type rho_;
};

} // namespace samples
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SAMPLES_DENSITY_MODE_HPP */
