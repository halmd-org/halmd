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
local lj_wrapper = {
    host = {
        [2] = halmd_wrapper.mdsim.host.forces.lj_2_
      , [3] = halmd_wrapper.mdsim.host.forces.lj_3_
      , potential = halmd_wrapper.mdsim.host.forces.lj_potential
    }
}
if halmd_wrapper.mdsim.gpu then
    lj_wrapper.gpu = {
        [2] = halmd_wrapper.mdsim.gpu.forces.lennard_jones_2_
      , [3] = halmd_wrapper.mdsim.gpu.forces.lennard_jones_3_
      , potential = halmd_wrapper.mdsim.gpu.forces.lj_potential
    }
end

local mdsim = {
    core = require("halmd.mdsim.core")
}
local device = require("halmd.device")
local args = require("halmd.options")
local assert = assert

module("halmd.mdsim.forces.lj", halmd.modules.register)

options = lj_wrapper.host.potential.options

--
-- construct lj module
--
function new()
    local dimension = assert(args.dimension)
    local cutoff = assert(args.cutoff)
    local epsilon = assert(args.epsilon)
    local sigma = assert(args.sigma)

    -- dependency injection
    local core = mdsim.core()
    local particle = assert(core.particle)
    local box = assert(core.box)

    if not device() then
        local potential = lj_wrapper.host.potential(particle.ntype, cutoff, epsilon, sigma)
        return lj_wrapper.host[dimension](potential, particle, box)
    end
    local potential = lj_wrapper.gpu.potential(particle.ntype, cutoff, epsilon, sigma)
    return lj_wrapper.gpu[dimension](potential, particle, box)
end
