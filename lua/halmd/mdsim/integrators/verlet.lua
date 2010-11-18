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
local verlet_wrapper = {
    host = {
        [2] = halmd_wrapper.mdsim.host.integrators.verlet_2_
      , [3] = halmd_wrapper.mdsim.host.integrators.verlet_3_
    }
}
if halmd_wrapper.mdsim.gpu then
    verlet_wrapper.gpu = {
        [2] = halmd_wrapper.mdsim.gpu.integrators.verlet_2_
      , [3] = halmd_wrapper.mdsim.gpu.integrators.verlet_3_
    }
end
local mdsim = {
  core = require("halmd.mdsim.core")
}
local po = halmd_wrapper.po
local device = require("halmd.device")
local assert = assert

module("halmd.mdsim.integrators.verlet", halmd.modules.register)

--
-- construct verlet module
--
function new(args)
    local timestep = args.timestep or 0.001 -- default value

    -- dependency injection
    local core = mdsim.core()
    local dimension = assert(core.dimension)
    local particle = assert(core.particle)
    local box = assert(core.box)

    if not device() then
        return verlet_wrapper.host[dimension](particle, box, timestep)
    end
    return verlet_wrapper.gpu[dimension](particle, box, timestep)
end

--
-- returns module description
--
function name()
    return "Velocity-Verlet integrator"
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    desc:add("timestep", po.float(), "integration timestep")
end
