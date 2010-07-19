/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <boost/mpl/if.hpp>

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/force_kernel.cuh>
#include <halmd/mdsim/gpu/forces/lj_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/gpu/blas/dsfloat.cuh>
#include <halmd/numeric/gpu/blas/symmetric.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/gpu/variant.cuh>

using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::numeric::gpu::blas;
using namespace halmd::utility::gpu;

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{
namespace lj_kernel
{

/** number of placeholders per neighbour list */
__constant__ unsigned int neighbour_size_;
/** neighbour list stride */
__constant__ unsigned int neighbour_stride_;
/** Lennard-Jones potential parameters */
texture<float4> ljparam_;
/** positions, types */
texture<float4> r_;
/** cuboid box edgle length */
__constant__ variant<map<pair<int_<3>, float3>, pair<int_<2>, float2> > > box_length_;

/**
 * Compute Lennard-Jones forces
 */
template <typename vector_type, typename gpu_vector_type>
__global__ void compute(
  gpu_vector_type* g_f,
  unsigned int* g_neighbour,
  float* g_en_pot,
  gpu_vector_type* g_virial)
{
    enum { dimension = vector_type::static_size };
    typedef typename vector_type::value_type value_type;
    unsigned int i = GTID;

    // load particle associated with this thread
    unsigned int type1;
    vector_type r1;
    tie(r1, type1) = untagged<vector_type>(tex1Dfetch(r_, i));

    // potential energy contribution
    float en_pot_ = 0;
    // virial contribution
    vector<float, (dimension - 1) * dimension / 2 + 1> virial_ = 0;
#ifdef USE_FORCE_DSFUN
    // force sum
    vector<dsfloat, dimension> f = 0;
#else
    vector_type f = 0;
#endif

    for (unsigned int k = 0; k < neighbour_size_; ++k) {
        // coalesced read from neighbour list
        unsigned int j = g_neighbour[k * neighbour_stride_ + i];
        // skip placeholder particles
        if (j == particle_kernel::PLACEHOLDER) {
            break;
        }

        // load particle
        unsigned int type2;
        vector_type r2;
        tie(r2, type2) = untagged<vector_type>(tex1Dfetch(r_, j));
        // Lennard-Jones potential parameters
        vector<float, 4> lj = tex1Dfetch(ljparam_, symmetric_matrix::lower_index(type1, type2));

        // particle distance vector
        vector_type r = r1 - r2;
        // enforce periodic boundary conditions
        box_kernel::reduce_periodic(r, static_cast<vector_type>(get<dimension>(box_length_)));
        // squared particle distance
        value_type rr = inner_prod(r, r);
        // enforce cutoff length
        if (rr >= lj[RR_CUT]) {
            continue;
        }

        // compute Lennard-Jones force in reduced units
        value_type rri = lj[SIGMA2] / rr;
        value_type ri6 = rri * rri * rri;
        value_type fval = 48 * lj[EPSILON] * rri * ri6 * (ri6 - 0.5f) / lj[SIGMA2];
        value_type en_pot = 4 * lj[EPSILON] * ri6 * (ri6 - 1) - lj[EN_CUT];

        // virial equation sum
        virial_ += 0.5f * fval * force_kernel::virial_tensor(rr, r);
        // potential energy contribution of this particle
        en_pot_ += 0.5f * en_pot;
        // force from other particle acting on this particle
        f += fval * r;
    }

    g_f[i] = static_cast<vector_type>(f);
    g_en_pot[i] = en_pot_;
    g_virial[i] = virial_;
}

} // namespace lj_kernel

template <int dimension>
lj_wrapper<dimension> const lj_wrapper<dimension>::kernel = {
    lj_kernel::r_
  , get<dimension>(lj_kernel::box_length_)
  , lj_kernel::neighbour_size_
  , lj_kernel::neighbour_stride_
  , lj_kernel::ljparam_
  , lj_kernel::compute<vector<float, dimension> >
};

template class lj_wrapper<3>;
template class lj_wrapper<2>;

}}} // namespace mdsim::gpu::forces

} // namespace halmd
