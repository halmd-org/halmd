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
local mdsim = {
    core = require("halmd.mdsim.core")
}
local neighbour_wrapper = {
    host = {
        [2] = halmd_wrapper.mdsim.host.neighbour_2_
      , [3] = halmd_wrapper.mdsim.host.neighbour_3_
    }
}
if halmd_wrapper.mdsim.gpu then
    neighbour_wrapper.gpu = {
        [2] = halmd_wrapper.mdsim.gpu.neighbour_2_
      , [3] = halmd_wrapper.mdsim.gpu.neighbour_3_
    }
end
local device = require("halmd.device")
local args = require("halmd.options")
local assert = assert

module("halmd.mdsim.neighbour", halmd.modules.register)

--
-- construct neighbour module
--
function new()
    -- dependency injection
    local core = mdsim.core()
    local particle = assert(core.particle)
    local box = assert(core.box)
    local force = assert(core.force)

    -- command line options
    local dimension = assert(args.dimension)
    local skin = assert(args.skin)

    if not device() then
        return neighbour_wrapper.host[dimension](particle, box, force, skin)
    end
    local cell_occupancy = assert(args.cell_occupancy)
    return neighbour_wrapper.gpu[dimension](particle, box, force, skin, cell_occupancy)
end

options = function(desc)
    neighbour_wrapper.host[2].options(desc)
    if neighbour_wrapper.gpu then
	neighbour_wrapper.gpu[2].options(desc)
    end
end
