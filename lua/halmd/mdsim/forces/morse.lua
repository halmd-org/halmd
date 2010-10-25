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
local morse_wrapper = {
    host = {
        [2] = halmd_wrapper.mdsim.host.forces.morse_2_
      , [3] = halmd_wrapper.mdsim.host.forces.morse_3_
      , potential = halmd_wrapper.mdsim.host.forces.morse_potential
    }
}
if halmd_wrapper.mdsim.gpu then
    morse_wrapper.gpu = {
        [2] = halmd_wrapper.mdsim.gpu.forces.morse_2_
      , [3] = halmd_wrapper.mdsim.gpu.forces.morse_3_
      , potential = halmd_wrapper.mdsim.gpu.forces.morse_potential
    }
end

local mdsim = {
    core = require("halmd.mdsim.core")
}
local device = require("halmd.device")
local args = require("halmd.options")
local assert = assert

module("halmd.mdsim.forces.morse", halmd.modules.register)

options = morse_wrapper.host.potential.options

--
-- construct module for Morse potential
--
function new()
    local dimension = assert(args.dimension)
    local cutoff = assert(args.cutoff)
    local epsilon = assert(args.epsilon)
    local sigma = assert(args.sigma)
    local rmin = assert(args.morse_minimum)

    -- dependency injection
    local core = mdsim.core()
    local particle = assert(core.particle)
    local box = assert(core.box)

    if not device() then
        local potential = morse_wrapper.host.potential(particle.ntype, cutoff, epsilon, sigma, rmin)
        return morse_wrapper.host[dimension](potential, particle, box)
    end
    local potential = morse_wrapper.gpu.potential(particle.ntype, cutoff, epsilon, sigma, rmin)
    return morse_wrapper.gpu[dimension](potential, particle, box)
end
