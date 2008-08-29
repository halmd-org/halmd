/* Global definitions
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

#ifndef MDSIM_CONFIG_HPP
#define MDSIM_CONFIG_HPP

#include "vector2d.hpp"
#include "vector3d.hpp"

/** positional coordinate dimensions */
#ifdef DIM_3D
static const unsigned int dimension = 3;
#else
static const unsigned int dimension = 2;
#endif

/** host vector type */
#ifdef DIM_3D
typedef vector3d<double> hvector;
#else
typedef vector2d<double> hvector;
#endif

#endif /* ! MDSIM_CONFIG_HPP */
