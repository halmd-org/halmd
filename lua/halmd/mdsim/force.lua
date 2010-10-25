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
local force_wrapper = {
    [2] = halmd_wrapper.mdsim.force_2_
  , [3] = halmd_wrapper.mdsim.force_3_
}
local forces = {
    lennard_jones = require("halmd.mdsim.forces.lennard_jones")
  , morse = require("halmd.mdsim.forces.morse")
  , power_law = require("halmd.mdsim.forces.power_law")
}
local args = require("halmd.options")
local assert = assert

module("halmd.mdsim.force", halmd.modules.register)

options = force_wrapper[2].options

--
-- construct force module
--
function new()
    local force = assert(args.force)
    return forces[force]()
end
