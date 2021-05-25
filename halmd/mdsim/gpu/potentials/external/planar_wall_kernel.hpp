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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_PLANAR_WALL_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_PLANAR_WALL_KERNEL_HPP

#include <halmd/utility/tuple.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {
namespace planar_wall_kernel {

/**
 * indices of potential parameters
 */
enum {
    EPSILON     /* interaction strength */
  , SIGMA       /* interaction range */
  , WETTING     /* wetting parameter */
  , CUTOFF      /* cutoff distance */
};

/**
 * planar_wall external potential
 */
template <int dimension>
class planar_wall
{
public:
    typedef fixed_vector<float, dimension> vector_type;

    /**
     * Construct planar_wall potential.
     */
    HALMD_GPU_ENABLED planar_wall(cudaTextureObject_t t_param_geometry, cudaTextureObject_t t_param_potential)
      : t_param_geometry_(t_param_geometry), t_param_potential_(t_param_potential)
    {}

    HALMD_GPU_ENABLED void fetch_param(unsigned int species);

     /**
     * Compute force and potential energy due to planar_wall walls.
     * Form of the potential is given here:
     * u_i(d)=epsilon_i*[(2/15)*(sigma_i/d)**9-wetting_i*(sigma_i/d)**3].
     */
    HALMD_GPU_ENABLED tuple<vector_type, float> operator()(vector_type const& r) const;


private:
    /** species of interacting particle */
    unsigned int species_;
    cudaTextureObject_t t_param_geometry_;
    cudaTextureObject_t t_param_potential_;
};

template <int dimension>
HALMD_GPU_ENABLED void planar_wall<dimension>::fetch_param(unsigned int species)
{
    species_ = species;
}

} // namespace planar_wall_kernel

struct planar_wall_wrapper
{
    static cuda::symbol<float> smoothing;
    static cuda::symbol<int> nwall;
};

} // namespace external
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_PLANAR_WALL_KERNEL_HPP */
