/* Lennard-Jones fluid simulation using CUDA
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef LJGPU_MDSIM_VARIANT_HPP
#define LJGPU_MDSIM_VARIANT_HPP

namespace ljgpu
{

enum mixture_type {
    // homogenous fluid
    UNARY,
    // binary mixture
    BINARY,
};

enum potential_type {
    // no smoothing
    C0POT,
    // C² potential smoothing
    C2POT,
};

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_VARIANT_HPP */
