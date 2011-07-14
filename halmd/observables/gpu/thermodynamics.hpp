/*
 * Copyright © 2010-2011  Felix Höfling and Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_THERMODYNAMICS_HPP
#define HALMD_OBSERVABLES_GPU_THERMODYNAMICS_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/observables/thermodynamics.hpp>

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
class thermodynamics
    : public observables::thermodynamics<dimension>
{
public:
    typedef observables::thermodynamics<dimension> _Base;
    typedef typename _Base::vector_type vector_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef typename _Base::box_type box_type;
    typedef typename _Base::clock_type clock_type;
    typedef mdsim::gpu::force<dimension, float_type> force_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<force_type> force;

    static void luaopen(lua_State* L);

    thermodynamics(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , boost::shared_ptr<clock_type> clock
      , boost::shared_ptr<force_type> force
    );

    virtual void prepare();
    virtual void sample(uint64_t step);

    virtual double en_kin();
    virtual vector_type v_cm();
    virtual double en_pot();
    virtual double virial();
    virtual double hypervirial();

private:
    /** functors for reduce kernels */
    algorithm::gpu::reduce<
        algorithm::gpu::sum_                    // reduce_transform
      , fixed_vector<float, dimension>          // input_type
      , float4                                  // coalesced_input_type
      , dsfloat                                 // output_type
      , dsfloat                                 // coalesced_output_type
      , double                                  // host_output_type
      , algorithm::gpu::square_                 // input_transform
    > sum_velocity_square_;

    algorithm::gpu::reduce<
        algorithm::gpu::sum_                    // reduce_transform
      , fixed_vector<float, dimension>          // input_type
      , float4                                  // coalesced_input_type
      , fixed_vector<dsfloat, dimension>        // output_type
      , fixed_vector<dsfloat, dimension>        // coalesced_output_type
      , vector_type                             // host_output_type
    > sum_velocity_vector_;

    algorithm::gpu::reduce<
        algorithm::gpu::sum_                    // reduce_transform
      , float                                   // input_type
      , float                                   // coalesced_input_type
      , dsfloat                                 // output_type
      , dsfloat                                 // coalesced_output_type
      , double                                  // host_output_type
    > sum_scalar_;

    typedef typename force_type::stress_tensor_type stress_tensor_type;
    typedef typename force_type::gpu_stress_tensor_type gpu_stress_tensor_type;
    algorithm::gpu::reduce<
        algorithm::gpu::sum_                    // reduce_transform
      , stress_tensor_type                      // input_type
      , gpu_stress_tensor_type                  // coalesced_input_type
      , dsfloat                                 // output_type
      , dsfloat                                 // coalesced_output_type
      , double                                  // host_output_type
      , algorithm::gpu::at_0                    // input_transform
    > sum_stress_tensor_diagonal_;
};

} // namespace observables
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_THERMODYNAMICS_HPP */
