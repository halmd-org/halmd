/*
 * Copyright © 2014 Sutapa Roy
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
#include <halmd/mdsim/gpu/potentials/external/slit_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {
namespace slit_kernel {

/** array of potential parameters for all particle species */
static texture<float4> param_potential_;
/** array of geometry parameters for all walls */
static texture<float4> param_geometry_;

/**
 * slit external potential
 */
template <int dimension>
class slit
{
public:
    typedef fixed_vector<float, dimension> vector_type;

    /**
     * Construct slit potential.
     *
     * Fetch parameters from texture cache for this particle species.
     */
    HALMD_GPU_ENABLED slit(unsigned int species) // FIXME , unsigned int nwall)
      : species_(species), nwall_(2)
    {}

     /**
     * Compute force and potential energy due to slit walls.
     * Form of the potential is given here:
     * u_i(d)=epsilon_i*[(2/15)*(sigma_i/d)**9-wetting_i*(sigma_i/d)**3].
     */
    HALMD_GPU_ENABLED tuple<vector_type, float> operator()(vector_type const& r) const
    {
        float en_pot = 0;
        vector_type force = 0;

        // loop over walls
        for (unsigned int i = 0; i < nwall_; ++i) {

            // fetch geometry parameters from texture cache
            vector_type surface_normal;
            float offset;
            tie(surface_normal, offset) <<= tex1Dfetch(param_geometry_, i);

            // compute absolute distance to wall i
            float d = inner_prod(r, surface_normal) - offset;
            
            // fetch potential parameters from texture cache
            float4 param_potential = tex1Dfetch(param_potential_, species_ * nwall_ + i);
            float epsilon = param_potential.x;
            float sigma = param_potential.y;
            float c = param_potential.z;

            // energy and force due to wall i
            float d3i = halmd::pow(sigma / d, 3);
            float d6i = float(2) / 15 * d3i * d3i;
            float eps_d3i = (d > 0 ? 1 : -1) * epsilon * d3i;
            float fval = 3 * eps_d3i * (3 * d6i - c) / d;
            force += fval * surface_normal;
            en_pot += eps_d3i * (d6i - c);
        }

        return make_tuple(force, en_pot);
    }


private:
    /** species of interacting particle */
    unsigned int species_;
    /** number of walls */
    unsigned int nwall_;

};

} // namespace slit_kernel

cuda::texture<float4> slit_wrapper::param_potential = slit_kernel::param_potential_;
cuda::texture<float4> slit_wrapper::param_geometry = slit_kernel::param_geometry_;

} // namespace external
} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using namespace potentials::external::slit_kernel;

template class external_wrapper<3, slit<3> >;
template class external_wrapper<2, slit<2> >;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
