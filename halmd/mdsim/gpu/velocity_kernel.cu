/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/gpu/velocity_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocity_kernel {

/** number of particles in simulation box */
static __constant__ unsigned int nbox_;

/**
 * Rescale magnitude of all velocities by 'factor'
 */
template <typename vector_type>
__global__ void rescale(
    float4* g_v
  , unsigned int size
  , dsfloat factor
)
{
    for (unsigned int i = GTID; i < nbox_; i += GTDIM) {
        vector_type v;
        unsigned int tag;
#ifdef USE_VERLET_DSFUN
        tie(v, tag) <<= tie(g_v[i], g_v[i + size]);
#else
        tie(v, tag) <<= g_v[i];
#endif
        v *= factor;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + size]) <<= tie(v, tag);
#else
        g_v[i] <<= tie(v, tag);
#endif
    }
}

/**
 * Shift all velocities by 'delta'
 */
template <typename vector_type>
__global__ void shift(
    float4* g_v
  , unsigned int size
  , fixed_vector<dsfloat, vector_type::static_size> delta
)
{
    for (unsigned int i = GTID; i < nbox_; i += GTDIM) {
        vector_type v;
        unsigned int tag;
#ifdef USE_VERLET_DSFUN
        tie(v, tag) <<= tie(g_v[i], g_v[i + size]);
#else
        tie(v, tag) <<= g_v[i];
#endif
        v += delta;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + size]) <<= tie(v, tag);
#else
        g_v[i] <<= tie(v, tag);
#endif
    }
}

/**
 * First shift, then rescale all velocities
 */
template <typename vector_type>
__global__ void shift_rescale(
    float4* g_v
  , unsigned int size
  , fixed_vector<dsfloat, vector_type::static_size> delta
  , dsfloat factor
)
{
    for (unsigned int i = GTID; i < nbox_; i += GTDIM) {
        vector_type v;
        unsigned int tag;
#ifdef USE_VERLET_DSFUN
        tie(v, tag) <<= tie(g_v[i], g_v[i + size]);
#else
        tie(v, tag) <<= g_v[i];
#endif
        v += delta;
        v *= factor;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + size]) <<= tie(v, tag);
#else
        g_v[i] <<= tie(v, tag);
#endif
    }
}

} // namespace velocity_kernel

template <int dimension>
velocity_wrapper<dimension> const velocity_wrapper<dimension>::kernel = {
#ifdef USE_VERLET_DSFUN
    velocity_kernel::rescale<fixed_vector<dsfloat, dimension> >
  , velocity_kernel::shift<fixed_vector<dsfloat, dimension> >
  , velocity_kernel::shift_rescale<fixed_vector<dsfloat, dimension> >
#else
    velocity_kernel::rescale<fixed_vector<float, dimension> >
  , velocity_kernel::shift<fixed_vector<float, dimension> >
  , velocity_kernel::shift_rescale<fixed_vector<float, dimension> >
#endif
  , velocity_kernel::nbox_
};

template class velocity_wrapper<3>;
template class velocity_wrapper<2>;

} // namespace mdsim
} // namespace gpu
} // namespace halmd
