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

#ifndef HALMD_MDSIM_FORCE_FLAGS_HPP
#define HALMD_MDSIM_FORCE_FLAGS_HPP

namespace halmd
{
namespace mdsim { namespace force_flags
{

// define compute flags for force modules
enum {
    potential_energy = 1 << 0
  , stress_tensor    = 1 << 1
  , hypervirial      = 1 << 2
};

}} // namespace mdsim::force_flags

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCE_FLAGS_HPP */
