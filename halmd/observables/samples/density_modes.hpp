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

#ifndef HALMD_OBSERVABLES_SAMPLES_DENSITY_MODES_HPP
#define HALMD_OBSERVABLES_SAMPLES_DENSITY_MODES_HPP

#include <boost/shared_ptr.hpp>
#include <complex>
#include <lua.hpp>
#include <vector>

#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd
{
namespace observables { namespace samples
{

/**
 *  data structure for storing Fourier modes of the particle density
 *
 *  @f$ \rho_{\vec q} = \sum_{i=1}^N \exp(\textrm{i}\vec q \cdot \vec r_i) @f$
 *  for each particle species
 */
template <int dimension>
class density_modes
{
public:
    typedef fixed_vector<double, dimension> vector_type;
    typedef std::complex<double> mode_type;

    /** density modes of a single particle species, one entry per wavevector */
    typedef std::vector<mode_type> mode_vector_type;
    /** list of density modes for all particle species */
    typedef std::vector<boost::shared_ptr<mode_vector_type> > mode_vector_vector_type;

    /**
     *  nested list of density modes, @f$ \text{rho[i][j]} = \rho_{\vec q}^{(i)} = @f$
     *  for wavevector @f$ \vec q = wavevector[j] @f$ and particle species @f$ i @f$
     */
    mode_vector_vector_type rho;
    /** simulation time when sample was taken */
    double time;

    static void luaopen(lua_State* L);

    density_modes(size_t n_species, size_t n_wavevectors);
};

}} // namespace observables::samples

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SAMPLES_DENSITY_MODES_HPP */
