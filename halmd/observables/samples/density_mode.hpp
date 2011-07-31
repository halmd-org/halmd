/*
 * Copyright © 2011  Felix Höfling
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

#include <boost/shared_ptr.hpp>
#include <complex>
#include <limits>
#include <lua.hpp>
#include <vector>

#include <halmd/mdsim/clock.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/raw_allocator.hpp>

namespace halmd {
namespace observables {
namespace samples {

/**
 *  data structure for storing Fourier modes of the particle density
 *
 *  @f$ \rho_{\vec q} = \sum_{i=1}^N \exp(\textrm{i}\vec q \cdot \vec r_i) @f$
 *  for each particle type
 */
template <int dimension>
class density_mode
{
private:
    typedef mdsim::clock clock_type;

public:
    typedef fixed_vector<double, dimension> vector_type;
    typedef std::complex<double> mode_type;
    typedef typename clock_type::step_type step_type;

    /** density modes of a single particle type, one entry per wavevector */
    typedef std::vector<mode_type, raw_allocator<mode_type> > mode_vector_type;
    /** list of density modes for all particle types */
    typedef std::vector<boost::shared_ptr<mode_vector_type> > mode_vector_vector_type;

    /**
     *  nested list of density modes, @f$ rho[i][j] = \rho_{\vec q}^{(i)} = @f$
     *  for wavevector @f$ \vec q = wavevector[j] @f$ and particle types @f$ i @f$
     */
    mode_vector_vector_type rho;
    /** simulation step when sample was taken */
    step_type step;

    static void luaopen(lua_State* L);
    static const char* class_name();

    /**
     * construct sample of given size
     *
     * @param ntype number of particle types
     * @param nq    total number of wavevectors
     */
    density_mode(unsigned int ntype, unsigned int nq);

    /**
     * Free shared pointers and re-allocate memory
     * if containers are shared with some other object.
     *
     * Values are not initialised.
     *
     * @param force if true then enforce reallocation
     */
    void reset(bool force=false);
};

template <int dimension>
inline density_mode<dimension>::density_mode(unsigned int ntype, unsigned int nq)
  // allocate sample pointers
  : rho(ntype)
  // initialise attributes
  , step(std::numeric_limits<step_type>::max())
{
    // allocate memory for each particle type
    for (unsigned int i = 0; i < ntype; ++i) {
        rho[i].reset(new mode_vector_type(nq));
    }
}

template <int dimension>
inline void density_mode<dimension>::reset(bool force)
{
    // free shared pointers and re-allocate memory
    for (size_t i = 0; i < rho.size(); ++i) {
        if (force || !rho[i].unique()) {
            rho[i].reset(new mode_vector_type(rho[i]->size()));
        }
    }
    // make time stamp invalid
    step = std::numeric_limits<step_type>::max();
}

} // namespace samples
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SAMPLES_DENSITY_MODE_HPP */
