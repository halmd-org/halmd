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
        [2] = halmd_wrapper.observables.host.samples.trajectory_2_double_
      , [3] = halmd_wrapper.observables.host.samples.trajectory_3_double_
    }
}
if halmd_wrapper.observables.gpu then
    trajectory_wrapper.gpu = {
        [2] = halmd_wrapper.observables.host.samples.trajectory_2_float_
      , [3] = halmd_wrapper.observables.host.samples.trajectory_3_float_
    }
end
local mdsim = {
    core = require("halmd.mdsim.core")
}
local device = require("halmd.device")
local assert = assert

module("halmd.observables.samples.trajectory", halmd.modules.register)

-- singleton
local sample

--
-- construct trajectory module
--
function new()
    if not sample then
        local core = assert(mdsim.core())
        local particle = assert(core.particle)
        local dimension = assert(core.dimension)

        local trajectory
        if device() then
            trajectory = assert(trajectory_wrapper.gpu[dimension])
        else
            trajectory = assert(trajectory_wrapper.host[dimension])
        end
        assert(trajectory)
        sample = trajectory(particle.ntypes)
    end
    return sample
end
