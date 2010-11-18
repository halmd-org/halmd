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
local mdsim = {
  core = require("halmd.mdsim.core")
}
local device = require("halmd.device")
local random = require("halmd.random")
local po = halmd_wrapper.po
local assert = assert

module("halmd.mdsim.integrators.verlet_nvt_andersen", halmd.modules.register)

--
-- construct verlet_nvt_andersen module
--
function new(args)
    local timestep = args.timestep or 0.01 -- default value
    local temperature = args.temperature or 1.12 -- default value
    local collision_rate = args.collision_rate or 10 -- default value

    -- dependency injection
    local core = mdsim.core()
    local dimension = assert(core.dimension)
    local particle = assert(core.particle)
    local box = assert(core.box)
    local random = assert(random())

    local andersen
    if device() then
        andersen = verlet_nvt_andersen_wrapper.gpu[dimension]
    else
        andersen = verlet_nvt_andersen_wrapper.host[dimension]
    end
    return andersen(particle, box, random, timestep, temperature, collision_rate)
end

--
-- returns module description
--
function name()
    return "Velocity-Verlet integrator with Andersen thermostat"
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    desc:add("timestep", po.float(), "integration timestep")
    desc:add("temperature", po.float(), "thermostat temperature")
    desc:add("collision-rate", po.float(), "collision rate for Andersen thermostat")
end
