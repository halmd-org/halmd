/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_DYNAMICS_MEAN_SQUARE_DISPLACEMENT_HPP
#define HALMD_OBSERVABLES_DYNAMICS_MEAN_SQUARE_DISPLACEMENT_HPP

#include <halmd/config.hpp>

namespace halmd {
namespace observables {
namespace dynamics {

template <typename vector_type>
struct mean_square_displacement
{
    typedef typename vector_type::value_type value_type;

    HALMD_GPU_ENABLED value_type operator()(vector_type const& r1, vector_type const& r2) const
    {
        // displacement of particle: R(t2) - R(t1)
        vector_type dr = r2 - r1;
        // return square displacement
        return inner_prod(dr, dr);
    }
};

} // namespace dynamics
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_DYNAMICS_MEAN_SQUARE_DISPLACEMENT_HPP */
