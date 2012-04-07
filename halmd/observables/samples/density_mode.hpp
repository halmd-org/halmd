/*
 * Copyright © 2011-2012  Felix Höfling and Peter Colberg
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
#include <limits>
#include <lua.hpp>
#include <stdexcept> // std::logic_error
#include <vector>

#include <halmd/mdsim/clock.hpp>
#include <halmd/utility/raw_allocator.hpp>

namespace halmd {
namespace observables {
namespace samples {

/**
 * Data structure for storing Fourier modes of the particle density.
 *
 * @f$ \rho_{\vec q} = \sum_{i=1}^N \exp(\textrm{i}\vec q \cdot \vec r_i) @f$
 */
template <int dimension>
class density_mode
{
public:
    typedef std::vector<std::complex<double>, raw_allocator<std::complex<double> > > mode_array_type;
    typedef typename mdsim::clock::step_type step_type;

    /**
     * Construct sample of given size.
     *
     * @param nq total number of wavevectors
     * @param step simulation step when sample is taken (optional)
     */
    density_mode(std::size_t nq, step_type step = std::numeric_limits<step_type>::max());

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
     * Returns simulation step when the sample was taken.
     */
    step_type step() const
    {
        if (step_ == std::numeric_limits<step_type>::max()) {
            throw std::logic_error("step not set in density modes sample");
        }
        return step_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** density modes */
    mode_array_type rho_;
    /** simulation step when sample was taken */
    step_type step_;
};

template <int dimension>
inline density_mode<dimension>::density_mode(std::size_t nq, step_type step)
  : rho_(nq)
  , step_(step)
{
}

} // namespace samples
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SAMPLES_DENSITY_MODE_HPP */
