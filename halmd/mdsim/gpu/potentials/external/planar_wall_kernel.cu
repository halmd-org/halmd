/*
 * Copyright © 2014-2015 Sutapa Roy
 * Copyright © 2014-2015 Felix Höfling
 * Copyright © 2020      Jaslo Ziska
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

#include <halmd/mdsim/gpu/forces/external_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/external/planar_wall_kernel.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function
#include <halmd/utility/tuple.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {
namespace planar_wall_kernel {

/** parameter for C2-smooth trunctation */
static __constant__ float smoothing_;
/** number of walls */
static __constant__ int nwall_;

template <int dimension>
__device__ tuple<typename planar_wall<dimension>::vector_type, float>
planar_wall<dimension>::operator()(vector_type const& r) const
{
    float en_pot = 0;
    vector_type force = 0;

    // loop over walls
    for (unsigned int i = 0; i < nwall_; ++i) {

        // fetch geometry parameters from texture cache
        vector_type surface_normal;
        float offset;
        tie(surface_normal, offset) <<= tex1Dfetch<float4>(t_param_geometry_, i);

        // compute absolute distance to wall i
        float d = inner_prod(r, surface_normal) - offset;

        // fetch potential parameters from texture cache
        fixed_vector<float, 4> param_potential = tex1Dfetch<float4>(t_param_potential_, species_ * nwall_ + i);
        float epsilon = param_potential[EPSILON];
        float sigma = param_potential[SIGMA];
        float w = param_potential[WETTING];
        float cutoff = param_potential[CUTOFF];

        // truncate interaction
        if (fabs(d) >= cutoff)
            continue;

        // cutoff energy due to wall i
        float dc3i = halmd::pow(sigma / cutoff, 3);
        float dc6i = float(2) / 15 * dc3i * dc3i;
        float eps_dc3i = epsilon * dc3i;
        float en_cut = eps_dc3i * (dc6i - w);

        // energy and force due to wall i
        float d3i = halmd::pow(sigma / d, 3);
        float d6i = float(2) / 15 * d3i * d3i;
        float eps_d3i = (d > 0 ? 1 : -1) * epsilon * d3i;
        float en_sub = eps_d3i * (d6i - w) - en_cut;
        float fval = 3 * eps_d3i * (3 * d6i - w) / d;

        // apply smooth truncation
        float dd = (fabs(d) - cutoff) / smoothing_;
        float x2 = dd * dd;
        float x4 = x2 * x2;
        float x4i = 1 / (1 + x4);
        float h0_r = x4 * x4i;

        // first derivative
        float h1_r = 4 * dd * x2 * x4i * x4i;

        // apply smoothing function to obtain C¹ force function
        fval = h0_r * fval - h1_r * en_sub / smoothing_;

        // apply smoothing function to obtain C² potential function
        en_sub = h0_r * en_sub;

        // accumulate force and potential energy
        force += fval * surface_normal;
        en_pot += en_sub;
    }
    return make_tuple(force, en_pot);
}

} // namespace planar_wall_kernel

cuda::symbol<float> planar_wall_wrapper::smoothing = planar_wall_kernel::smoothing_;
cuda::symbol<int> planar_wall_wrapper::nwall = planar_wall_kernel::nwall_;

} // namespace external
} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using namespace potentials::external::planar_wall_kernel;

template class external_wrapper<3, planar_wall<3> >;
template class external_wrapper<2, planar_wall<2> >;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
