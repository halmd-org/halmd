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

#ifndef HALMD_MDSIM_FORCE_KERNEL_HPP
#define HALMD_MDSIM_FORCE_KERNEL_HPP

#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {

/**
 * Trace and off-diagonal elements of distance tensor
 */
template <typename float_type>
HALMD_GPU_ENABLED typename type_traits<3, float_type>::stress_tensor_type
make_stress_tensor(float_type rr, fixed_vector<float_type, 3> const& r)
{
    typename type_traits<3, float_type>::stress_tensor_type v;
    v[0] = rr;
    v[1] = r[1] * r[2];
    v[2] = r[2] * r[0];
    v[3] = r[0] * r[1];
    return v;
}

template <typename float_type>
HALMD_GPU_ENABLED typename type_traits<2, float_type>::stress_tensor_type
make_stress_tensor(float_type rr, fixed_vector<float_type, 2> const& r)
{
    typename type_traits<2, float_type>::stress_tensor_type v;
    v[0] = rr;
    v[1] = r[0] * r[1];
    return v;
}

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_FORCE_KERNEL_HPP */
