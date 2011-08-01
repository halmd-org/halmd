/*
 * Copyright © 2008-2011  Peter Colberg
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

#include <halmd/observables/gpu/dynamics/mean_square_displacement_kernel.hpp>
#include <halmd/observables/gpu/dynamics/tagged_particle.cuh>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

template <typename vector_type>
struct mean_square_displacement_kernel
{
    typedef typename vector_type::value_type value_type;

    HALMD_GPU_ENABLED value_type operator()(vector_type const& r1, vector_type const& r2) const
    {
        vector_type dr = r2 - r1;
        return inner_prod(dr, dr);
    }
};

template <int dimension, unsigned int threads>
mean_square_displacement_wrapper<dimension, threads> const
mean_square_displacement_wrapper<dimension, threads>::wrapper = {
    accumulate<mean_square_displacement_kernel<vector_type>, threads>
};

// explicit template instantiation
template class mean_square_displacement_wrapper<3, 512>;
template class mean_square_displacement_wrapper<2, 512>;
template class mean_square_displacement_wrapper<3, 256>;
template class mean_square_displacement_wrapper<2, 256>;
template class mean_square_displacement_wrapper<3, 128>;
template class mean_square_displacement_wrapper<2, 128>;
template class mean_square_displacement_wrapper<3, 64>;
template class mean_square_displacement_wrapper<2, 64>;
template class mean_square_displacement_wrapper<3, 32>;
template class mean_square_displacement_wrapper<2, 32>;

} // namespace dynamics
} // namespace gpu
} // namespace observables
} // namespace halmd
