/*
 * Copyright © 2014 Felix Höfling
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

#include <halmd/mdsim/gpu/forces/external_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/external/harmonic_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {
namespace harmonic_kernel {

/** array of potential parameters for all particle species */
static texture<float4> param_;

/**
 * harmonic external potential
 */
template <int dimension>
class harmonic
{
public:
    typedef fixed_vector<float, dimension> vector_type;

    /**
     * Construct harmonic potential.
     *
     * Fetch parameters from texture cache for this particle species.
     */
    HALMD_GPU_ENABLED harmonic(unsigned int species)
    {
        tie(offset_, stiffness_) <<= tex1Dfetch(param_, species);
    }

    /**
     * Compute force and potential for interaction.
     *
     * @param r   particle position reduced to periodic box
     * @returns   tuple of force vector @f$ -\nabla U(\vec r) @f$ and potential @f$ U(\vec r) @f$
     */
    HALMD_GPU_ENABLED tuple<vector_type, float> operator()(vector_type const& r) const
    {
        vector_type dr = r - offset_;
        vector_type force = stiffness_ * dr;
        float en_pot = inner_prod(force, dr) / 2;

        return make_tuple(force, en_pot);
    }

private:
    /** potential parameters for given particle species */
    vector_type offset_;
    float stiffness_;
};

} // namespace harmonic_kernel

cuda::texture<float4> harmonic_wrapper::param = harmonic_kernel::param_;

} // namespace external
} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using namespace potentials::external::harmonic_kernel;

template class external_wrapper<3, harmonic<3> >;
template class external_wrapper<2, harmonic<2> >;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
