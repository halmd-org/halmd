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
local po = halmd_wrapper.po
local assert = assert

module("halmd.mdsim.forces.morse", halmd.modules.register)

--
-- construct Morse module
--
function new(args)
    local cutoff = args.cutoff or { 2.5, 2.5, 2.5 } -- default value
    local epsilon = args.epsilon or { 1.0, 1.5, 0.5 } -- default value
    local sigma = args.sigma or { 1.0, 0.8, 0.88 } -- default value
    local minimum = args.morse_minimum or { 2.0, 1.6, 1.78 } -- default value

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

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    -- FIXME desc:add("cutoff", po.float_array(), "truncate potential at cutoff radius")
    -- FIXME desc:add("epsilon", po.float_array(), "potential well depths")
    -- FIXME desc:add("sigma", po.float_array(), "collision diameters")
    desc:add("morse-minimum", po.float_array(), "positions of potential minimum for interactions")
    -- FIXME desc:add("smooth", po.float_array(), "C²-potential smoothing factor")
end
