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
template <int dimension, typename float_type, typename ptr_type>
__global__ void rescale(
    ptr_type g_v
  , unsigned int nparticle
  , unsigned int size
  , dsfloat factor
)
{
    for (unsigned int i = GTID; i < nparticle; i += GTDIM) {
        fixed_vector<float_type, dimension> v;
        unsigned int id;
        tie(v, id) <<= g_v[i];
        v *= factor;
        g_v[i] <<= tie(v, id);
    }
}

/**
 * Rescale magnitude of velocities of group by 'factor'
 */
template <int dimension, typename float_type, typename ptr_type>
__global__ void rescale_group(
    ptr_type g_v
  , unsigned int const* g_group
  , unsigned int nparticle
  , unsigned int size
  , dsfloat factor
)
{
    for (unsigned int n = GTID; n < nparticle; n += GTDIM) {
        unsigned int i = g_group[n];
        fixed_vector<float_type, dimension> v;
        unsigned int id;
        tie(v, id) <<= g_v[i];
        v *= factor;
        g_v[i] <<= tie(v, id);
    }
}

/**
 * Shift all velocities by 'delta'
 */
template <int dimension, typename float_type, typename ptr_type>
__global__ void shift(
    ptr_type g_v
  , unsigned int nparticle
  , unsigned int size
  , fixed_vector<dsfloat, dimension> delta
)
{
    for (unsigned int i = GTID; i < nparticle; i += GTDIM) {
        fixed_vector<float_type, dimension> v;
        unsigned int id;
        tie(v, id) <<= g_v[i];
        v += delta;
        g_v[i] <<= tie(v, id);
    }
}

/**
 * Shift velocities of group by 'delta'
 */
template <int dimension, typename float_type, typename ptr_type>
__global__ void shift_group(
    ptr_type g_v
  , unsigned int const* g_group
  , unsigned int nparticle
  , unsigned int size
  , fixed_vector<dsfloat, dimension> delta
)
{
    for (unsigned int n = GTID; n < nparticle; n += GTDIM) {
        unsigned int i = g_group[n];
        fixed_vector<float_type, dimension> v;
        unsigned int id;
        tie(v, id) <<= g_v[i];
        v += delta;
        g_v[i] <<= tie(v, id);
    }
}

/**
 * First shift, then rescale all velocities
 */
template <int dimension, typename float_type, typename ptr_type>
__global__ void shift_rescale(
    ptr_type g_v
  , unsigned int nparticle
  , unsigned int size
  , fixed_vector<dsfloat, dimension> delta
  , dsfloat factor
)
{
    for (unsigned int i = GTID; i < nparticle; i += GTDIM) {
        fixed_vector<float_type, dimension> v;
        unsigned int id;
        tie(v, id) <<= g_v[i];
        v += delta;
        v *= factor;
        g_v[i] <<= tie(v, id);
    }
}

/**
 * First shift, then rescale velocities of group
 */
template <int dimension, typename float_type, typename ptr_type>
__global__ void shift_rescale_group(
    ptr_type g_v
  , unsigned int const* g_group
  , unsigned int nparticle
  , unsigned int size
  , fixed_vector<dsfloat, dimension> delta
  , dsfloat factor
)
{
    for (unsigned int n = GTID; n < nparticle; n += GTDIM) {
        unsigned int i = g_group[n];
        fixed_vector<float_type, dimension> v;
        unsigned int id;
        tie(v, id) <<= g_v[i];
        v += delta;
        v *= factor;
        g_v[i] <<= tie(v, id);
    }
}

} // namespace velocity_kernel

template <int dimension, typename float_type>
velocity_wrapper<dimension, float_type> velocity_wrapper<dimension, float_type>::kernel = {
    velocity_kernel::rescale<dimension, float_type, ptr_type>
  , velocity_kernel::rescale_group<dimension, float_type, ptr_type>
  , velocity_kernel::shift<dimension, float_type, ptr_type>
  , velocity_kernel::shift_group<dimension, float_type, ptr_type>
  , velocity_kernel::shift_rescale<dimension, float_type, ptr_type>
  , velocity_kernel::shift_rescale_group<dimension, float_type, ptr_type>
};

template class velocity_wrapper<3, float>;
template class velocity_wrapper<2, float>;
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class velocity_wrapper<3, dsfloat>;
template class velocity_wrapper<2, dsfloat>;
#endif

} // namespace mdsim
} // namespace gpu
} // namespace halmd
