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
local sampler_wrapper = {
    [2] = halmd_wrapper.sampler_2_
  , [3] = halmd_wrapper.sampler_3_
}
local mdsim = {
    core = require("halmd.mdsim.core")
}
local args = require("halmd.options")
local assert = assert
local math = math

module("halmd.sampler", halmd.modules.register)

options = sampler_wrapper[2].options

local sampler -- singleton instance

--
-- construct sampler module
--
function new()
    -- dependency injection
    local core = mdsim.core()
    local integrator = assert(core.integrator)

    -- command line options
    local dimension = assert(args.dimension)
    local sampling_state_vars = assert(args.sampling_state_vars)
    local sampling_trajectory = assert(args.sampling_trajectory)
    local steps = assert(args.steps)
    local time = args.time -- optional

    if not sampler then
        local wrapper = sampler_wrapper[dimension]
        if time then
            steps = math.floor((time / integrator.timestep) + 0.5)
        end
        sampler = wrapper(core, steps, sampling_state_vars, sampling_trajectory)
    end
    return sampler
end
