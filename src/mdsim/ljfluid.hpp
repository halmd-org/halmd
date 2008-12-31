/* Lennard-Jones fluid simulation
 *
 * Copyright (C) 2008  Peter Colberg
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

#ifndef MDSIM_LJFLUID_HPP
#define MDSIM_LJFLUID_HPP

#ifdef USE_CUDA
# if defined(USE_NEIGHBOUR) || !defined(USE_CELL)
#  include "ljfluid_gpu_nbr.hpp"
# else
#  include "ljfluid_gpu_cell.hpp"
# endif
#else
# include "ljfluid_host.hpp"
#endif

#endif /* ! MDSIM_LJFLUID_HPP */
