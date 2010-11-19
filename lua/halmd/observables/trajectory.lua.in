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
local trajectory_wrapper = {
    host = {
        [2] = halmd_wrapper.observables.host.trajectory_2_
      , [3] = halmd_wrapper.observables.host.trajectory_3_
    }
}
if halmd_wrapper.observables.gpu then
    trajectory_wrapper.gpu = {
        [2] = halmd_wrapper.observables.gpu.host.trajectory_2_
      , [3] = halmd_wrapper.observables.gpu.host.trajectory_3_
    }
end
local mdsim = {
    core = require("halmd.mdsim.core")
}
local observables = {
    samples = {
        trajectory = require("halmd.observables.samples.trajectory")
    }
}
local device = require("halmd.device")
local assert = assert

module("halmd.observables.trajectory", halmd.modules.register)

--
-- construct trajectory module
--
function new(args)
    -- dependency injection
    local core = mdsim.core()
    local dimension = assert(core.dimension)
    local particle = assert(core.particle)
    local box = assert(core.box)
    local sample = assert(observables.samples.trajectory())

    local trajectory
    if device() then
        trajectory = assert(trajectory_wrapper.gpu[dimension])
    else
        trajectory = assert(trajectory_wrapper.host[dimension])
    end
    return trajectory(sample, particle, box)
end
