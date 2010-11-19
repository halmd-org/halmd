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
local h5 = halmd_wrapper.h5
local po = halmd_wrapper.po
local assert = assert
local math = math

module("halmd.sampler", halmd.modules.register)

local sampler -- singleton instance

--
-- construct sampler module
--
function new(args)
    -- dependency injection
    local core = mdsim.core()
    local dimension = assert(core.dimension)
    local integrator = assert(core.integrator)

    -- command line options
    local state_vars = args.state_vars or 25 -- default value
    local trajectory = args.trajectory or 0 -- default value
    local steps = args.steps or 10000 -- default value
    local time = args.time -- optional

    if not sampler then
        local wrapper = sampler_wrapper[dimension]
        if time then
            steps = math.floor((time / integrator.timestep) + 0.5)
        end
        sampler = wrapper(core, steps, state_vars, trajectory)
    end
    return sampler
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc, globals)
    globals:add("steps", po.uint64(), "number of simulation steps")
    globals:add("time", po.float():conflicts("steps"), "total simulation time")
    desc:add("state-vars", po.uint(), "sample interval for macroscopic state variables")
    desc:add("trajectory", po.uint(), "sample interval for trajectory")
end

--
-- write module parameters to HDF5 group
--
-- @param sampler module instance
-- @param group HDF5 group
--
function write_parameters(sampler, group, globals)
    globals:write_attribute("steps", h5.uint64(), sampler.steps)
    globals:write_attribute("time", h5.float(), sampler.total_time)
    group:write_attribute("state_vars", h5.uint(), sampler.statevars_interval)
    group:write_attribute("trajectory", h5.uint(), sampler.trajectory_interval)
end
