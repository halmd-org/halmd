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

#ifndef HALMD_OBSERVABLES_UTILITY_WAVEVECTORS_HPP
#define HALMD_OBSERVABLES_UTILITY_WAVEVECTORS_HPP

#include <map>
#include <vector>

#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd
{
namespace observables { namespace utility
{

/**
 * construct set of wavevector shells compatible with the reciprocal
 * lattice and a list of wavenumbers
 */

template <int dimension>
class wavevectors
{
public:
    typedef fixed_vector<double, dimension> vector_type;

    wavevectors(
        std::vector<double> const& wavenumbers
      , vector_type const& box_length
      , double tolerance
      , unsigned int max_count
    );

    //! returns tolerance on wavevector magnitude
    double tolerance() const
    {
        return tolerance_;
    }

    //! returns maximum count of wavevectors per wavenumber
    unsigned int maximum_count() const
    {
        return max_count_;
    }

    //! returns list of wavevectors
    std::multimap<double, vector_type> const& values() const
    {
        return wavevectors_;
    }

    //! returns wavenumber grid
    std::vector<double> const& wavenumbers() const
    {
        return wavenumbers_;
    }

protected:
    /** wavenumber grid */
    std::vector<double> wavenumbers_;
    /** tolerance of wavevector magnitudes (relative error) */
    double tolerance_;
    /** maximum number of wavevectors per wavenumber */
    double max_count_;
    /**
     * list of wavevectors grouped by their magnitude,
     * the keys equal wavenumbers_ (or a subset of)
     */
    std::multimap<double, vector_type> wavevectors_;
};

}} // namespace observables::utility

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_UTILITY_WAVEVECTORS_HPP */
