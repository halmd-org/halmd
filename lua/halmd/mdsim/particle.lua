--
-- Copyright © 2010  Peter Colberg
--
-- This file is part of HALMD.
--
-- HALMD is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with this program.  If not, see <http://www.gnu.org/licenses/>.
--

require("halmd.modules")

-- grab environment
local particle_wrapper = {
    host = {
        [2] = halmd_wrapper.mdsim.host.particle_2_
      , [3] = halmd_wrapper.mdsim.host.particle_3_
    }
  , [2] = halmd_wrapper.mdsim.particle_2_
  , [3] = halmd_wrapper.mdsim.particle_3_
}
if halmd_wrapper.mdsim.gpu then
    particle_wrapper.gpu = {
        [2] = halmd_wrapper.mdsim.gpu.particle_2_
      , [3] = halmd_wrapper.mdsim.gpu.particle_3_
    }
end
local device = require("halmd.device")
local args = require("halmd.options")
local assert = assert

module("halmd.mdsim.particle", halmd.modules.register)

-- option/parameter namespace
namespace = "box"

options = particle_wrapper[2].options

--
-- construct particle module
--
function new()
    local dimension = assert(args.dimension)
    local npart = assert(args.particles)

    if not device() then
        return particle_wrapper.host[dimension](npart)
    end
    return particle_wrapper.gpu[dimension](device(), npart)
end
