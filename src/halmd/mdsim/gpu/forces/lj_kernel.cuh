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

#ifndef HALMD_MDSIM_GPU_FORCES_LJ_KERNEL_CUH
#define HALMD_MDSIM_GPU_FORCES_LJ_KERNEL_CUH

namespace halmd { namespace mdsim { namespace gpu { namespace forces { namespace lj_kernel
{

//
// Lennard Jones potential parameter indices
//
enum {
    /** potential well depths in MD units */
    EPSILON,
    /** square of cutoff length */
    RR_CUT,
    /** square of pair separation */
    SIGMA2,
    /** potential energy at cutoff length in MD units */
    EN_CUT,
};

}}}}} // namespace halmd::mdsim::gpu::forces::lj_kernel

#endif /* ! HALMD_MDSIM_GPU_FORCES_LJ_KERNEL_CUH */
