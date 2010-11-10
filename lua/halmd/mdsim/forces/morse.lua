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
    host = halmd_wrapper.mdsim.host.forces.morse
}
if halmd_wrapper.mdsim.gpu then
    morse_wrapper.gpu = halmd_wrapper.mdsim.gpu.forces.morse
end

local mdsim = {
    core = require("halmd.mdsim.core")
}
local device = require("halmd.device")
local args = require("halmd.options")
local assert = assert
local hooks = require("halmd.hooks")

module("halmd.mdsim.forces.morse", halmd.modules.register)

options = morse_wrapper.host.options

--
-- construct Morse module
--
function new()
    local cutoff = assert(args.cutoff)
    local epsilon = assert(args.epsilon)
    local sigma = assert(args.sigma)
    local minimum = assert(args.morse_minimum)

    local core = mdsim.core()
    local particle = assert(core.particle)

    local morse
    if device() then
        morse = assert(morse_wrapper.gpu)
    else
        morse = assert(morse_wrapper.host)
    end
    return morse(particle.ntype, cutoff, epsilon, sigma, minimum)
end
