/*
 * Copyright © 2011  Felix Höfling
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

#ifndef HALMD_OBSERVABLES_GPU_FIELDS_DENSITY_KERNEL_HPP
#define HALMD_OBSERVABLES_GPU_FIELDS_DENSITY_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace fields {

struct density_wrapper
{
    /** compute number density from cell list and cell volume */
    cuda::function<void (unsigned int const*, double*, float, unsigned int)> compute;

    static density_wrapper kernel;
};

} // namespace fields
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_FIELDS_DENSITY_KERNEL_HPP */
