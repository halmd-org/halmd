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
    host = {
        [2] = halmd_wrapper.mdsim.host.sampler.trajectory_2_
      , [3] = halmd_wrapper.mdsim.host.sampler.trajectory_3_
    }
}
if halmd_wrapper.mdsim.gpu then
    sampler_wrapper.gpu = {
        [2] = halmd_wrapper.mdsim.gpu.sampler.trajectory_host_2_
      , [3] = halmd_wrapper.mdsim.gpu.sampler.trajectory_host_3_
    }
end
local mdsim = {
    core = require("halmd.mdsim.core")
}
local args = require("halmd.options")
local assert = assert

module("halmd.mdsim.samples.trajectory", halmd.modules.register)

--
-- construct trajectory module
--
function new()
    -- dependency injection
    local core = mdsim.core()
    local particle = assert(core.particle)
    local box = assert(core.box)

    -- command line options
    local dimension = assert(args.dimension)
    local backend = assert(args.backend)

    return sampler_wrapper[backend][dimension](particle, box, core)
end
