--
-- Copyright © 2010  Peter Colberg and Felix Höfling
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
local mdsim = {
  core = require("halmd.mdsim.core")
}
local verlet_nvt_andersen_wrapper = {
    host = {
        [2] = halmd_wrapper.mdsim.host.integrators.verlet_nvt_andersen_2_
      , [3] = halmd_wrapper.mdsim.host.integrators.verlet_nvt_andersen_3_
    }
}
if halmd_wrapper.mdsim.gpu then
    verlet_nvt_andersen_wrapper.gpu = {
        [2] = halmd_wrapper.mdsim.gpu.integrators.verlet_nvt_andersen_2_
      , [3] = halmd_wrapper.mdsim.gpu.integrators.verlet_nvt_andersen_3_
    }
end
local device = require("halmd.device")
local random = {
    gpu = require("halmd.gpu.random")
  , host = require("halmd.host.random")
}
local args = require("halmd.options")
local assert = assert

module("halmd.mdsim.integrators.verlet_nvt_andersen", halmd.modules.register)

options = verlet_nvt_andersen_wrapper.host[2].options

--
-- construct verlet_nvt_andersen module
--
function new()
    local dimension = assert(args.dimension)
    local timestep = assert(args.timestep)
    local temperature = assert(args.temperature)
    local collision_rate = args.andersen_collision_rate or 10

    -- dependency injection
    local core = mdsim.core()
    local particle = assert(core.particle)
    local box = assert(core.box)

    if not device() then
        return verlet_nvt_andersen_wrapper.host[dimension](
            particle, box, random.host(), timestep, temperature, collision_rate
        )
    end
    return verlet_nvt_andersen_wrapper.gpu[dimension](
        particle, box, random.gpu(), timestep, temperature, collision_rate
    )
end
