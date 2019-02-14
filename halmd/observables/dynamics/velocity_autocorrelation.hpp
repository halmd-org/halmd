/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_OBSERVABLES_DYNAMICS_VELOCITY_AUTOCORRELATION_HPP
#define HALMD_OBSERVABLES_DYNAMICS_VELOCITY_AUTOCORRELATION_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace observables {
namespace dynamics {

template <unsigned int dimension, typename float_type>
struct velocity_autocorrelation
{
    typedef fixed_vector<float_type, dimension> vector_type;

    HALMD_GPU_ENABLED float_type operator()(vector_type const& v1, vector_type const& v2) const
    {
        return inner_prod(v1, v2);
    }
};

} // namespace dynamics
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_DYNAMICS_VELOCITY_AUTOCORRELATION_HPP */
