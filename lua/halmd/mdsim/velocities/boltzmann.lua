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
local boltzmann_wrapper = {
    host = {
        [2] = halmd_wrapper.mdsim.host.velocities.boltzmann_2_
      , [3] = halmd_wrapper.mdsim.host.velocities.boltzmann_3_
    }
}
if halmd_wrapper.mdsim.gpu then
    boltzmann_wrapper.gpu = {
        [2] = halmd_wrapper.mdsim.gpu.velocities.boltzmann_2_
      , [3] = halmd_wrapper.mdsim.gpu.velocities.boltzmann_3_
    }
end
local mdsim = {
    core = require("halmd.mdsim.core")
}
local random = {
    gpu = require("halmd.gpu.random")
  , host = require("halmd.host.random")
}
local device = require("halmd.device")
local po = halmd_wrapper.po
local assert = assert

module("halmd.mdsim.velocities.boltzmann", halmd.modules.register)

--
-- construct boltzmann module
--
function new(args)
    local temperature = args.temperature or 1.12 -- default value

    -- dependency injection
    local core = mdsim.core()
    local dimension = assert(core.dimension)
    local particle = assert(core.particle)

    if not device() then
        return boltzmann_wrapper.host[dimension](particle, random.host(), temperature)
    end
    return boltzmann_wrapper.gpu[dimension](particle, random.gpu(), temperature)
end

--
-- returns module description
--
function name()
    return "Boltzmann distribution"
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    desc:add("temperature", po.float(), "Boltzmann distribution temperature")
end
