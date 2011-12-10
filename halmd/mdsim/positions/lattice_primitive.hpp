/*
 * Copyright © 2008  Peter Colberg
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

#ifndef HALMD_MDSIM_POSITIONS_LATTICE_PRIMITIVE_HPP
#define HALMD_MDSIM_POSITIONS_LATTICE_PRIMITIVE_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {

/**
 * place particles on a face centered cubic lattice (fcc)
 */
struct fcc_lattice_primitive
{
    template <typename float_type>
    HALMD_GPU_ENABLED void operator()(
        fixed_vector<float_type, 3>& r
      , fixed_vector<unsigned int, 3> const& nsite
      , unsigned int site
    ) const
    {
        r[0] = ((site >> 2) % nsite[0]) + ((site ^ (site >> 1)) & 1) / float_type(2);
        r[1] = ((site >> 2) / nsite[0] % nsite[1]) + (site & 1) / float_type(2);
        r[2] = ((site >> 2) / nsite[0] / nsite[1]) + (site & 2) / float_type(4);
    }

    template <typename float_type>
    HALMD_GPU_ENABLED void operator()(
        fixed_vector<float_type, 2>& r
      , fixed_vector<unsigned int, 2> const& nsite
      , unsigned int site
    ) const
    {
        r[0] = ((site >> 1) % nsite[0]) + (site & 1) / float_type(2);
        r[1] = ((site >> 1) / nsite[0]) + (site & 1) / float_type(2);
    }
};

/**
 * place particles on a simple cubic lattice (sc)
 */
struct sc_lattice_primitive
{
    template <typename float_type>
    HALMD_GPU_ENABLED void operator()(
        fixed_vector<float_type, 3>& r
      , fixed_vector<unsigned int, 3> const& nsite
      , unsigned int site
    ) const
    {
        r[0] = (site % nsite[0]) + float_type(0.5);
        r[1] = (site / nsite[0] % nsite[1]) + float_type(0.5);
        r[2] = (site / nsite[0] / nsite[1]) + float_type(0.5);
    }

    template <typename float_type>
    HALMD_GPU_ENABLED void operator()(
        fixed_vector<float_type, 2>& r
      , fixed_vector<unsigned int, 2> const& nsite
      , unsigned int site
    ) const
    {
        r[0] = (site % nsite[0]) + float_type(0.5);
        r[1] = (site / nsite[0]) + float_type(0.5);
    }
};

} // namespace halmd

#endif /* ! HALMD_MDSIM_POSITIONS_LATTICE_PRIMITIVE_HPP */
