/*
 * Copyright © 2016       Manuel Dibak
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

#ifndef HALMD_MDSIM_GPU_FORCES_ASYMMETRIC_TRUNC_KERNEL_CUH
#define HALMD_MDSIM_GPU_FORCES_ASYMMETRIC_TRUNC_KERNEL_CUH

#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/forces/asymmetric_trunc_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {
namespace asymmetric_trunc_kernel {

/** positions, types */
static texture<float4> r2_;
/** orientation, nothing */
static texture<float4> u2_;

/**
 * Compute pair forces, potential energy, and stress tensor for all particles
 */
template <
    bool do_aux               //< compute auxiliary variables in addition to force
  , typename vector_type
  , typename potential_type
  , typename gpu_vector_type
  , typename trunc_type
>
__global__ void compute(
    float4 const* g_r1
  , float4 const* g_u1
  , gpu_vector_type* g_f
  , gpu_vector_type* g_tau
  , unsigned int const* g_neighbour
  , unsigned int neighbour_size
  , unsigned int neighbour_stride
  , float* g_en_pot
  , float* g_stress_pot
  , unsigned int ntype1
  , unsigned int ntype2
  , vector_type box_length
  , trunc_type const trunc
  , bool force_zero
  , bool torque_zero    // FIXME flag is unneeded
  , float aux_weight
)
{
    enum { dimension = vector_type::static_size };
    typedef typename vector_type::value_type value_type;
    typedef typename type_traits<dimension, float>::pseudo_vector_type torque_type;
    typedef typename type_traits<dimension, float>::stress_tensor_type stress_tensor_type;
    unsigned int i = GTID;

    // load particle associated with this thread
    unsigned int type1, nothing;
    vector_type r1, u1;
    tie(r1, type1) <<= g_r1[i];
    tie(u1, nothing) <<= g_u1[i];

    // contribution to potential energy
    float en_pot_ = 0;
    // contribution to stress tensor
    stress_tensor_type stress_pot = 0;
#ifdef USE_FORCE_DSFUN
    // force sum
    fixed_vector<dsfloat, dimension> f_ = 0;
    torque_type tau_ = 0;
#else
    vector_type f_ = 0;
    torque_type tau_ = 0;
#endif
    for (unsigned int k = 0; k < neighbour_size; ++k) {

      // coalesced read from neighbour list
      unsigned int j = g_neighbour[k * neighbour_stride + i];
      // skip placeholder particles
      if (j == particle_kernel::placeholder) {
          break;
      }
      // load particle
      unsigned int type2, nothing;
      vector_type r2, u2;
      tie(r2, type2) <<= tex1Dfetch(r2_, j);
      tie(u2, nothing) <<= tex1Dfetch(u2_, j);
      // pair potential
      potential_type const potential(type1, type2, ntype1, ntype2);

      // particle distance vector
      vector_type r = r1 - r2;
      // enforce periodic boundary conditions
      box_kernel::reduce_periodic(r, box_length);
      // squared particle distance
      value_type rr = inner_prod(r, r);
      // enforce cutoff length
      if (!potential.within_range(r, u1, u2)) {
          continue;
      }
      value_type en_pot;
      vector_type f;
      torque_type tau;
      tie(f, tau, en_pot) = potential(r, u1, u2);
      f_ += f;
      tau_ += tau;
      // apply smoothing function to force and potential
      value_type fval = sqrtf( inner_prod(f, f) );
      //trunc(sqrt(rr), sqrt(potential.rr_cut()), fval, en_pot);

      if (do_aux) {
          // potential energy contribution of this particle
          en_pot_ += aux_weight * en_pot;
          // contribution to stress tensor from this particle
          stress_pot += aux_weight * fval * make_stress_tensor(r);
      }
    }

    // add old force and auxiliary variables if not zero
    if (!force_zero) {
        f_ += static_cast<vector_type>(g_f[i]);
        tau_ += static_cast<torque_type>(g_tau[i]);
        if (do_aux) {
            en_pot_ += g_en_pot[i];
            stress_pot += read_stress_tensor<stress_tensor_type>(g_stress_pot + i, GTDIM);
        }
    }

    // write results to global memory
    g_f[i] = static_cast<vector_type>(f_);
    g_tau[i] = static_cast<torque_type>(tau_);
    if (do_aux) {
        g_en_pot[i] = en_pot_;
        write_stress_tensor(g_stress_pot + i, stress_pot, GTDIM);
    }
}

} // namespace asymmetric_trunc_kernel

template <int dimension, typename potential_type, typename trunc_type>
asymmetric_trunc_wrapper<dimension, potential_type, trunc_type> const
asymmetric_trunc_wrapper<dimension, potential_type, trunc_type>::kernel = {
    asymmetric_trunc_kernel::compute<false, fixed_vector<float, dimension>, potential_type>
  , asymmetric_trunc_kernel::compute<true, fixed_vector<float, dimension>, potential_type>
  , asymmetric_trunc_kernel::r2_
  , asymmetric_trunc_kernel::u2_
};

} // namespace mdsim
} // namespace gpu
} // namespace forces
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_ASYMMETRIC_TRUNC_KERNEL_CUH */
