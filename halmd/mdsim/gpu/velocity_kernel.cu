/*
 * Copyright © 2008-2010, 2012 Peter Colberg
 * Copyright © 2010 Felix Höfling
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

#include <halmd/mdsim/gpu/velocity_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocity_kernel {

/**
 * Rescale magnitude of all velocities by 'factor'
 */
template <typename vector_type>
__global__ void rescale(
    float4* g_v
  , unsigned int nparticle
  , unsigned int size
  , dsfloat factor
)
{
    for (unsigned int i = GTID; i < nparticle; i += GTDIM) {
        vector_type v;
        unsigned int id;
#ifdef USE_VERLET_DSFUN
        tie(v, id) <<= tie(g_v[i], g_v[i + size]);
#else
        tie(v, id) <<= g_v[i];
#endif
        v *= factor;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + size]) <<= tie(v, id);
#else
        g_v[i] <<= tie(v, id);
#endif
    }
}

/**
 * Rescale magnitude of velocities of group by 'factor'
 */
template <typename vector_type>
__global__ void rescale_group(
    float4* g_v
  , unsigned int const* g_group
  , unsigned int nparticle
  , unsigned int size
  , dsfloat factor
)
{
    for (unsigned int n = GTID; n < nparticle; n += GTDIM) {
        unsigned int i = g_group[n];
        vector_type v;
        unsigned int id;
#ifdef USE_VERLET_DSFUN
        tie(v, id) <<= tie(g_v[i], g_v[i + size]);
#else
        tie(v, id) <<= g_v[i];
#endif
        v *= factor;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + size]) <<= tie(v, id);
#else
        g_v[i] <<= tie(v, id);
#endif
    }
}

/**
 * Shift all velocities by 'delta'
 */
template <typename vector_type>
__global__ void shift(
    float4* g_v
  , unsigned int nparticle
  , unsigned int size
  , fixed_vector<dsfloat, vector_type::static_size> delta
)
{
    for (unsigned int i = GTID; i < nparticle; i += GTDIM) {
        vector_type v;
        unsigned int id;
#ifdef USE_VERLET_DSFUN
        tie(v, id) <<= tie(g_v[i], g_v[i + size]);
#else
        tie(v, id) <<= g_v[i];
#endif
        v += delta;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + size]) <<= tie(v, id);
#else
        g_v[i] <<= tie(v, id);
#endif
    }
}

/**
 * Shift velocities of group by 'delta'
 */
template <typename vector_type>
__global__ void shift_group(
    float4* g_v
  , unsigned int const* g_group
  , unsigned int nparticle
  , unsigned int size
  , fixed_vector<dsfloat, vector_type::static_size> delta
)
{
    for (unsigned int n = GTID; n < nparticle; n += GTDIM) {
        unsigned int i = g_group[n];
        vector_type v;
        unsigned int id;
#ifdef USE_VERLET_DSFUN
        tie(v, id) <<= tie(g_v[i], g_v[i + size]);
#else
        tie(v, id) <<= g_v[i];
#endif
        v += delta;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + size]) <<= tie(v, id);
#else
        g_v[i] <<= tie(v, id);
#endif
    }
}

/**
 * First shift, then rescale all velocities
 */
template <typename vector_type>
__global__ void shift_rescale(
    float4* g_v
  , unsigned int nparticle
  , unsigned int size
  , fixed_vector<dsfloat, vector_type::static_size> delta
  , dsfloat factor
)
{
    for (unsigned int i = GTID; i < nparticle; i += GTDIM) {
        vector_type v;
        unsigned int id;
#ifdef USE_VERLET_DSFUN
        tie(v, id) <<= tie(g_v[i], g_v[i + size]);
#else
        tie(v, id) <<= g_v[i];
#endif
        v += delta;
        v *= factor;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + size]) <<= tie(v, id);
#else
        g_v[i] <<= tie(v, id);
#endif
    }
}

/**
 * First shift, then rescale velocities of group
 */
template <typename vector_type>
__global__ void shift_rescale_group(
    float4* g_v
  , unsigned int const* g_group
  , unsigned int nparticle
  , unsigned int size
  , fixed_vector<dsfloat, vector_type::static_size> delta
  , dsfloat factor
)
{
    for (unsigned int n = GTID; n < nparticle; n += GTDIM) {
        unsigned int i = g_group[n];
        vector_type v;
        unsigned int id;
#ifdef USE_VERLET_DSFUN
        tie(v, id) <<= tie(g_v[i], g_v[i + size]);
#else
        tie(v, id) <<= g_v[i];
#endif
        v += delta;
        v *= factor;
#ifdef USE_VERLET_DSFUN
        tie(g_v[i], g_v[i + size]) <<= tie(v, id);
#else
        g_v[i] <<= tie(v, id);
#endif
    }
}

} // namespace velocity_kernel

template <int dimension>
velocity_wrapper<dimension> const velocity_wrapper<dimension>::kernel = {
#ifdef USE_VERLET_DSFUN
    velocity_kernel::rescale<fixed_vector<dsfloat, dimension> >
  , velocity_kernel::rescale_group<fixed_vector<dsfloat, dimension> >
  , velocity_kernel::shift<fixed_vector<dsfloat, dimension> >
  , velocity_kernel::shift_group<fixed_vector<dsfloat, dimension> >
  , velocity_kernel::shift_rescale<fixed_vector<dsfloat, dimension> >
  , velocity_kernel::shift_rescale_group<fixed_vector<dsfloat, dimension> >
#else
    velocity_kernel::rescale<fixed_vector<float, dimension> >
  , velocity_kernel::rescale_group<fixed_vector<float, dimension> >
  , velocity_kernel::shift<fixed_vector<float, dimension> >
  , velocity_kernel::shift_group<fixed_vector<float, dimension> >
  , velocity_kernel::shift_rescale<fixed_vector<float, dimension> >
  , velocity_kernel::shift_rescale_group<fixed_vector<float, dimension> >
#endif
};

template class velocity_wrapper<3>;
template class velocity_wrapper<2>;

} // namespace mdsim
} // namespace gpu
} // namespace halmd
