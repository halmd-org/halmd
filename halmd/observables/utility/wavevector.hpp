/*
 * Copyright © 2011-2013 Felix Höfling
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

#include <lua.hpp>
#include <utility>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace observables {
namespace utility {

/**
 * construct set of wavevector shells compatible with the reciprocal
 * lattice of the periodic simulation box and a list of wavenumbers
 */

template <int dimension>
class wavevector
{
public:
    typedef fixed_vector<double, dimension> vector_type;
    typedef std::vector<double> wavenumber_array_type;
    typedef std::vector<vector_type> wavevector_array_type;
    typedef std::vector<std::pair<std::size_t, std::size_t>> shell_array_type;

    static void luaopen(lua_State* L);

    // constructor
    wavevector(
        std::vector<double> const& wavenumber
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
    unsigned int max_count() const
    {
        return max_count_;
    }

    //! returns list of wavevectors
    wavevector_array_type const& value() const
    {
        return wavevector_;
    }

    //! returns wavenumber grid
    wavenumber_array_type const& wavenumber() const
    {
        return wavenumber_;
    }

    /*
     * returns list of wavevector shells
     *
     * The values are index ranges on the wavevector array, the order of shells
     * coincides with the wavenumber grid.
     */
    shell_array_type const& shell() const
    {
        return shell_;
    }

protected:
    /** wavenumber grid */
    wavenumber_array_type wavenumber_;
    /** edge lengths of simulation box */
    vector_type box_length_;
    /** tolerance of wavevector magnitudes (relative error) */
    double tolerance_;
    /** maximum number of wavevectors per wavenumber */
    double max_count_;
    // list of wavevectors grouped by their magnitude in ascending order
    wavevector_array_type wavevector_;
    // list of wavevector shells
    shell_array_type shell_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace observables
} // namespace utility
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_UTILITY_WAVEVECTOR_HPP */
