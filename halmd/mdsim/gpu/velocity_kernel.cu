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
template <typename float_type, int dimension>
__global__ void rescale(
    typename type_traits<4, float_type>::gpu::ptr_type g_v
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
template <typename float_type, int dimension>
__global__ void rescale_group(
    typename type_traits<4, float_type>::gpu::ptr_type g_v
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
template <typename float_type, int dimension>
__global__ void shift(
    typename type_traits<4, float_type>::gpu::ptr_type g_v
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
template <typename float_type, int dimension>
__global__ void shift_group(
    typename type_traits<4, float_type>::gpu::ptr_type g_v
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
template <typename float_type, int dimension>
__global__ void shift_rescale(
    typename type_traits<4, float_type>::gpu::ptr_type g_v
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
template <typename float_type, int dimension>
__global__ void shift_rescale_group(
    typename type_traits<4, float_type>::gpu::ptr_type g_v
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

template <typename float_type, int dimension>
velocity_wrapper<float_type, dimension> const velocity_wrapper<float_type, dimension>::kernel = {
    velocity_kernel::rescale<float_type, dimension>
  , velocity_kernel::rescale_group<float_type, dimension>
  , velocity_kernel::shift<float_type, dimension>
  , velocity_kernel::shift_group<float_type, dimension>
  , velocity_kernel::shift_rescale<float_type, dimension>
  , velocity_kernel::shift_rescale_group<float_type, dimension>
};

template class velocity_wrapper<float, 3>;
template class velocity_wrapper<float, 2>;
template class velocity_wrapper<dsfloat, 3>;
template class velocity_wrapper<dsfloat, 2>;

} // namespace mdsim
} // namespace gpu
} // namespace halmd
