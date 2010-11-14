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
local thermodynamics_wrapper = {
    host = {
        [2] = halmd_wrapper.observables.host.thermodynamics_2_
      , [3] = halmd_wrapper.observables.host.thermodynamics_3_
    }
  , [2] = halmd_wrapper.observables.thermodynamics_2_
  , [3] = halmd_wrapper.observables.thermodynamics_3_
}
if halmd_wrapper.observables.gpu then
    thermodynamics_wrapper.gpu = {
        [2] = halmd_wrapper.observables.gpu.thermodynamics_2_
      , [3] = halmd_wrapper.observables.gpu.thermodynamics_3_
    }
end
local mdsim = {
    core = require("halmd.mdsim.core")
}
local device = require("halmd.device")
local assert = assert

module("halmd.observables.thermodynamics", halmd.modules.register)

--
-- construct thermodynamics module
--
function new(args)
    -- dependency injection
    local core = mdsim.core()
    local dimension = assert(core.dimension)
    local particle = assert(core.particle)
    local box = assert(core.box)
    local force = assert(core.force)

    if not device() then
        return thermodynamics_wrapper.host[dimension](particle, box, force)
    end
    return thermodynamics_wrapper.gpu[dimension](particle, box, force)
end
