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

#ifndef HALMD_OBSERVABLES_UTILITY_WAVEVECTOR_HPP
#define HALMD_OBSERVABLES_UTILITY_WAVEVECTOR_HPP

#include <utility>
#include <lua.hpp>
#include <vector>

#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd
{
namespace observables { namespace utility
{

/**
 * construct set of wavevector shells compatible with the reciprocal
 * lattice and a list of wavenumbers
 *
 * The result in value() is a sorted container of key/value pairs
 * (wavenumber, wavevector).
 */

template <int dimension>
class wavevector
{
public:
    typedef fixed_vector<double, dimension> vector_type;
    typedef std::vector<std::pair<double, vector_type> > map_type;

    static void luaopen(lua_State* L);

    // construct class with list of wavenumbers
    wavevector(
        std::vector<double> const& wavenumber
      , vector_type const& box_length
      , double tolerance
      , unsigned int max_count
    );

    // construct class with upper limit on wavenumber,
    // the grid is linearly spaced starting with the smallest value
    // that is compatible with the extents of the simulation box
    wavevector(
        double max_wavenumber
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
    map_type const& value() const
    {
        return wavevector_;
    }

    //! returns wavenumber grid
    std::vector<double> const& wavenumber() const
    {
        return wavenumber_;
    }

protected:
    // common part of constructors
    void init_();

    /** wavenumber grid */
    std::vector<double> wavenumber_;
    /** edge lengths of simulation box */
    vector_type box_length_;
    /** tolerance of wavevector magnitudes (relative error) */
    double tolerance_;
    /** maximum number of wavevectors per wavenumber */
    double max_count_;
    /**
     * list of wavevectors grouped by their magnitude,
     * the keys equal wavenumber_ (or a subset of)
     */
    map_type wavevector_;
};

}} // namespace observables::utility

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_UTILITY_WAVEVECTOR_HPP */
