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
local power_law_wrapper = {
    host = {
        [2] = halmd_wrapper.mdsim.host.forces.power_law_2_
      , [3] = halmd_wrapper.mdsim.host.forces.power_law_3_
    }
}
if halmd_wrapper.mdsim.gpu then
    power_law_wrapper.gpu = {
        [2] = halmd_wrapper.mdsim.gpu.forces.power_law_2_
      , [3] = halmd_wrapper.mdsim.gpu.forces.power_law_3_
    }
end
local mdsim = {
  core = require("halmd.mdsim.core")
}
local args = require("halmd.options")
local assert = assert

module("halmd.mdsim.forces.power_law", halmd.modules.register)

options = power_law_wrapper.host[2].options

--
-- construct power_law module
--
function new()
    -- dependency injection
    local core = mdsim.core()
    local particle = assert(core.particle)
    local box = assert(core.box)

    -- command line options
    local dimension = assert(args.dimension)
    local backend = assert(args.backend)
    local cutoff = assert(args.cutoff)
    local epsilon = assert(args.epsilon)
    local sigma = assert(args.sigma)
    local index = assert(args.power_law_index)

    local power_law = power_law_wrapper[backend][dimension]
    return power_law(particle, box, index, cutoff, epsilon, sigma)
end
