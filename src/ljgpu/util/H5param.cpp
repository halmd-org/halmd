/* HDF5 parameter group
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

#include <ljgpu/util/H5param.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <ljgpu/util/exception.hpp>
#include <ljgpu/version.h>
#include <stdint.h>

namespace ljgpu
{

/**
 * initialize HDF5 parameter group
 */
H5param::H5param(H5::Group param) : param(param)
{
    H5xx::group node(param.createGroup("program"));
    // write program info attributes
    node["name"] = PROGRAM_NAME;
    node["version"] = PROGRAM_VERSION;
    node["variant"] = PROGRAM_VARIANT;
}

} // namespace ljgpu
