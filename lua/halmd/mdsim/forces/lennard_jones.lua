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
local lennard_jones_wrapper = {
    host = halmd_wrapper.mdsim.host.forces.lennard_jones
}
if halmd_wrapper.mdsim.gpu then
    lennard_jones_wrapper.gpu = halmd_wrapper.mdsim.gpu.forces.lennard_jones
end

local mdsim = {
    core = require("halmd.mdsim.core")
}
local device = require("halmd.device")
local args = require("halmd.options")
local assert = assert
local hooks = require("halmd.hooks")

module("halmd.mdsim.forces.lennard_jones", halmd.modules.register)

options = lennard_jones_wrapper.host.options

--
-- construct Lennard-Jones module
--
function new()
    local cutoff = assert(args.cutoff)
    local epsilon = assert(args.epsilon)
    local sigma = assert(args.sigma)

    local core = mdsim.core()
    local particle = assert(core.particle)

    local lennard_jones
    if device() then
        lennard_jones = assert(lennard_jones_wrapper.gpu)
    else
        lennard_jones = assert(lennard_jones_wrapper.host)
    end
    return lennard_jones(particle.ntype, cutoff, epsilon, sigma)
end
