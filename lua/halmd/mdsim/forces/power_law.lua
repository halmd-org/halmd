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
local power_law_wrapper = {
    host = halmd_wrapper.mdsim.host.forces.power_law
}
if halmd_wrapper.mdsim.gpu then
    power_law_wrapper.gpu = halmd_wrapper.mdsim.gpu.forces.power_law
end

local mdsim = {
    core = require("halmd.mdsim.core")
}
local device = require("halmd.device")
local args = require("halmd.options")
local assert = assert
local hooks = require("halmd.hooks")

module("halmd.mdsim.forces.power_law", halmd.modules.register)

options = power_law_wrapper.host.options

--
-- construct power law module
--
function new()
    local index = assert(args.power_law_index)
    local cutoff = assert(args.cutoff)
    local epsilon = assert(args.epsilon)
    local sigma = assert(args.sigma)

    local core = mdsim.core()
    local particle = assert(core.particle)

    local power_law
    if device() then
        power_law = assert(power_law_wrapper.gpu)
    else
        power_law = assert(power_law_wrapper.host)
    end
    return power_law(particle.ntype, index, cutoff, epsilon, sigma)
end
